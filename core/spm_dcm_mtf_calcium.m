function [S, Hz] = spm_dcm_mtf_calcium(P, M, U)
%==========================================================================
% Transfer function for 2-state calcium DCM
% FORMAT [S, Hz] = spm_dcm_mtf_calcium(P, M, U)
%
% INPUT:
%   P - Parameters structure with A, C, etc.
%   M - Model structure with f, g, x, etc.
%   U - Input structure (optional)
%
% OUTPUT:
%   S  - Transfer function [nFreq x nRegions x nRegions]
%   Hz - Frequency vector
%
% THEORY:
%   S(ω, i→j) = [C_out * (jωI - J)^{-1} * B_in](ω)
%
%   where:
%     J       - Jacobian of full system (2n×2n)
%     B_in    - Input maps to neuronal states (2n×n)
%     C_out   - Observation Jacobian dgdx from M.g at x0 (n×2n)
%               M1 Hill:   C_out(:,xCa) = diag(alpha .* dF/dCa)
%               M2 linear: C_out(:,xCa) = diag(alpha)
%               Consistent with manuscript hat_kappa = kappa*tau*alpha*F'(c*)
%     ω       - Angular frequency (2πf)
%
% The calcium model has 2 states per region:
%   x(:,1) - Neuronal activity
%   x(:,2) - Calcium concentration (observable)
%
% CRITICAL DIFFERENCE from spm_dcm_mtf:
%   - Standard spm_dcm_mtf assumes 1 state per region
%   - This function handles 2 states per region correctly
%   - Computes Jacobian respecting [n×2] state structure
%
%--------------------------------------------------------------------------
%   Park & ChatGPT (2025)

n = size(P.A,1);
m = n*2;

% Frequency grid
if isfield(M,'Hz'), Hz = M.Hz(:);
else
    maxHz=5; %for calcium imaging
    dt = spm_field(M,'dt',1/2);
    N  = spm_field(M,'N',32);
    Hz = linspace(0.01,maxHz,N)';   % 1–5 Hz default
end
Nf = length(Hz);

% Steady state
if isfield(M,'x') && ~isempty(M.x)
    x0 = full(M.x);
else
    x0 = zeros(n,2);
end

f  = M.f; if ~isa(f,'function_handle'), f = str2func(f); end

% Operating point for Jacobian evaluation.
% M.x is set by spm_dcm_calcium_priors: x(:,1)=0, x(:,2)=x2_ss≈0.12.
% The v2 voltage mapping guarantees sigmaV(x1=0)=0.5 for any V0,
% so x1=0 is always the correct neuronal operating point.
%
% Only substitute when M.x was not provided at all (genuinely all-zero).
% Do NOT override when x0(:,1)=0 but x0(:,2) is already set by priors —
% that is the normal operating condition and Ca must be kept at its
% prior steady-state, not forced to an arbitrary value.
%
% Ca steady-state: solve f2(x1=0, x2_ss)=0 from the state function.
% f2 = kappa*sigmaV - x2/tau_Ca + beta_x is linear in x2, so
% two evaluations at x2=0 and x2=1 give the exact root via interpolation.
% This is parameter-general (works for any P, not just the prior).
if all(x0(:) == 0)
    x0(:,1) = 0.0;
    fa = reshape(feval(f, [zeros(n,1), zeros(n,1)], [], P, M), n, 2);
    fb = reshape(feval(f, [zeros(n,1),  ones(n,1)], [], P, M), n, 2);
    f2a = fa(:,2);  f2b = fb(:,2);
    den = f2b - f2a;
    den(abs(den) < 1e-12) = -1;            % guard: flat eqn → x2_ss=0
    x0(:,2) = max(-f2a ./ den, 0.01);     % x2_ss, clamped positive
end

% Get analytic Jacobian from state function
[f0, J_analytic] = f(x0,[],P,M);

% Use analytic Jacobian if available, otherwise compute numerically
if ~isempty(J_analytic) && isnumeric(J_analytic)
    if iscell(J_analytic)
        % Convert cell Jacobian to full matrix
        J = [J_analytic{1,1}, J_analytic{1,2}; ...
             J_analytic{2,1}, J_analytic{2,2}];
        J = full(J);
    else
        J = full(J_analytic);
    end
else
    % Numerical Jacobian as fallback
    xv = spm_vec(x0);  J = zeros(m,m);  dx = 1e-6;
    for j=1:m
        x1v = xv;
        x1v(j) = x1v(j) + dx;
        x1 = spm_unvec(x1v, x0);
        f1 = f(x1,[],P,M);
        J(:,j) = (spm_vec(f1)-spm_vec(f0))/dx;
    end
end

% I/O maps
B_in = sparse(m,n);
B_in(1:n,1:n) = speye(n);


% Observation Jacobian at operating point — parallels spm_dcm_mtf.m lines 83-87.
% dgdx from M.g: Hill gives [0|diag(alpha*dF/dCa)]; linear gives [0|diag(alpha)]
if isfield(M,'g') && ~isempty(M.g)
    g = M.g;
    if ~isa(g,'function_handle'), g = str2func(g); end
    if nargout(g) >= 2
        [~, dgdx] = feval(g, x0, [], P, M);
    else
        dgdx = spm_diff(g, x0, [], P, M, 1);
    end
    C_out = sparse(dgdx);
else
    C_out = sparse(n,m);
    for i=1:n, C_out(i,n+i)=1; end   % fallback: unit obs
end

% Transfer function
S = zeros(Nf,n,n);
for k=1:Nf
    omega = 2*pi*Hz(k);
    H = ((1i*omega*eye(m) - J) \ eye(m));
    S(k,:,:) = C_out * H * B_in;
end
end
