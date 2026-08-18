function [f, dfdx, D, dfdu] = spm_fx_calcium(x, u, P, M)
%==========================================================================
% spm_fx_calcium  —  State equations for DCM with calcium dynamics
%
% FORMAT [f,dfdx,D,dfdu] = spm_fx_calcium(x,u,P,M)
%
% States:
%   x(:,1) - neuronal activity
%   x(:,2) - intracellular calcium concentration
%
% Parameters (in P):
%   A,B,C,D    : effective connectivity parameters
%   tau_Ca     : log decay constant
%   kappa_Ca   : log influx scaling
%   beta_Ca    : log baseline influx
%   R          : log sigmoid slope
%   V0         : resting membrane voltage (mV); σ(x1=0)=0.5 by construction
%
% Options (in M):
%   M.use_numeric = true → compute Jacobian numerically
%
% Voltage mapping (v2, 2026-06-06):
%   V_eff = V0 + x1*(V_max - V0)
%   At x1=0: V_eff = V0  →  sigmaV = 0.5 exactly, for any V0.
%   This guarantees the DCM linearisation point (x1=0) sits at the
%   sigmoid's operating point regardless of V0.  Prior: pE.V0 = 0.
%
%   Previous mapping (v0, archived as spm_fx_calcium_v0.m):
%   V_eff = -70 + ((x1+1)/2)*120  (= -10 + 60*x1)
%   Required V0 prior to be manually set to -10 mV to achieve σ(x1=0)=0.5.
%
%--------------------------------------------------------------------------
% Park (2025)
%==========================================================================

%--------------------------------------------------------------------------
% Default handling
%--------------------------------------------------------------------------
if nargin < 4, M = struct(); end
P.A = full(P.A); P.B = full(P.B); P.C = full(P.C); P.D = full(P.D);
n  = size(P.A,1);
nu = numel(u(:));

% ensure state as [n×2]
if size(x,2) < 2
    if numel(x) == 2*n, x = reshape(full(x), n, 2);
    else, x = full(x); x(:,2) = 0; end
else
    x = full(x);
end
x1 = x(:,1); x2 = x(:,2);

%--------------------------------------------------------------------------
% Calcium parameters (global scaling)
%--------------------------------------------------------------------------
H       = [0.8 0.02 0.26 0.2];  % [tau_base, beta_base, kappa_base, R_base]
tau_Ca  = H(1) * exp(P.tau_Ca);
beta_x  = H(2) * exp(P.beta_Ca(:));
kappa   = H(3) * exp(P.kappa_Ca);
R       = H(4) * exp(P.R);

V0      = P.V0(:); if isempty(V0), V0 = zeros(n,1); end

%--------------------------------------------------------------------------
% Voltage mapping: x1 → V_eff  (v2: operating point guaranteed at x1=0)
%--------------------------------------------------------------------------
V_max  = 50;
V_eff  = V0 + x1 .* (V_max - V0);   % x1=0 → V_eff=V0 → sigmaV=0.5
sigmaV = 1 ./ (1 + exp(-R*(V_eff - V0)));

%--------------------------------------------------------------------------
% Effective connectivity EE = A + Σ u_i B_i + Σ x_i D_i
%--------------------------------------------------------------------------
EE = P.A;
for i = 1:size(P.B,3)
    if ~isempty(u) && i <= nu, EE = EE + u(i)*P.B(:,:,i); end
end
for i = 1:size(P.D,3)
    EE = EE + x1(i)*P.D(:,:,i);
end

% Self-inhibition correction (SPM convention)
SE = diag(EE);
EE = EE - diag(exp(SE)/2 + SE);

%--------------------------------------------------------------------------
% Flow equations
%--------------------------------------------------------------------------
f = zeros(n,2);

if ~isempty(P.C) && ~isempty(u)
    f(:,1) = EE * x1 + full(P.C) * u(:);
else
    f(:,1) = EE * x1;
end

f(:,2) = kappa .* sigmaV - (1/tau_Ca) .* x2 + beta_x;
f = f(:);

%--------------------------------------------------------------------------
% Return if only flows needed
%--------------------------------------------------------------------------
if nargout < 2
    D = 1;
    return
end

use_numeric = false;
if isfield(M,'use_numeric') && M.use_numeric
    use_numeric = true;
elseif isfield(P,'D') && ~isempty(P.D)
    if any(abs(spm_vec(P.D)) > 1e-8)
        use_numeric = true;
    end
end
%==========================================================================
% === ANALYTIC JACOBIAN ===================================================
%==========================================================================
if ~use_numeric

    dfdx = cell(2,2);
    dfdu = cell(2,1);

    %--------------------------------------------------------------
    % (1) Neuronal block dfdx{1,1}
    %--------------------------------------------------------------
    dfdx{1,1} = EE;
    for i = 1:size(P.D,3)
        Di = P.D(:,:,i);
        Di = Di + diag( (diag(EE) - 1) .* diag(Di) );
        dfdx{1,1}(:,i) = dfdx{1,1}(:,i) + Di * x1;
    end

    %--------------------------------------------------------------
    % (2) Neuronal → calcium coupling  (dV/dx1 = V_max - V0)
    %--------------------------------------------------------------
    sigp = sigmaV .* (1 - sigmaV) .* (R .* (V_max - V0));
    dfdx{2,1} = diag(kappa(:) .* sigp(:));

    %--------------------------------------------------------------
    % (3) Calcium self-decay
    %--------------------------------------------------------------
    dfdx{1,2} = sparse(n,n);
    dfdx{2,2} = -speye(n) * (1/tau_Ca);

    %--------------------------------------------------------------
    % (4) Concatenate Jacobian
    %--------------------------------------------------------------
    dfdx = spm_cat(dfdx);

    %--------------------------------------------------------------
    % (5) Input Jacobian
    %--------------------------------------------------------------
    dfdu{1,1} = P.C;
    for i = 1:size(P.B,3)
        Bi = P.B(:,:,i);
        Bi = Bi + diag( (diag(EE) - 1) .* diag(Bi) );
        dfdu{1,1}(:,i) = dfdu{1,1}(:,i) + Bi * x1;
    end
    dfdu{2,1} = sparse(n, numel(u(:)));
    dfdu = spm_cat(dfdu);

    %--------------------------------------------------------------
    % (6) Delay operator
    %--------------------------------------------------------------
    D = 1;

%==========================================================================
% === NUMERICAL JACOBIAN ==================================================
%==========================================================================
else
    f0   = f;
    xv   = spm_vec(x);
    nx   = numel(xv);
    du   = u(:);
    nu   = numel(du);
    dfdx = zeros(numel(f0), nx);
    h    = 1e-6;

    for k = 1:nx
        x1p = spm_unvec(xv, x); x1p(k) = x1p(k) + h;
        f1  = spm_fx_calcium(x1p, u, P, rmfield(M,'use_numeric'));
        dfdx(:,k) = (f1(:) - f0(:)) / h;
    end

    if nargout > 3
        dfdu = zeros(numel(f0), nu);
        for j = 1:nu
            uj = du; uj(j) = uj(j) + h;
            f1 = spm_fx_calcium(x, uj, P, rmfield(M,'use_numeric'));
            dfdu(:,j) = (f1(:) - f0(:)) / h;
        end
    else
        dfdu = sparse(size(f0,1),0);
    end

    dfdx = sparse(dfdx);
    dfdu = sparse(dfdu);
    D = 1;
end
end
