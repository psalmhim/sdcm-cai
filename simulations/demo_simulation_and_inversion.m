%==========================================================================
% demo_simulation_and_inversion.m
%
% Self-contained demo: simulates 4-node population calcium imaging data
% from a known ground-truth connectivity matrix, then inverts the
% simulated data twice with spm_dcm_calcium_csd.m — once under the
% default Hill (saturating) observation model, once under the linear
% observation model — and compares recovered connectivity and free
% energy between the two.
%
% To switch observation models, set:
%   DCM.options.custom_g = @spm_gx_calcium_linear;   % linear
% (omit, or leave empty, for the default Hill model)
%
% Requires SPM25 (or SPM12) and this repository's core/ folder on the
% MATLAB path.
%==========================================================================

clear; close all; rng(1);

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core')));

%--------------------------------------------------------------------------
% Ground-truth connectivity (4 nodes) — same structure used in the paper's
% Fig. 2 4-node recovery simulation
%--------------------------------------------------------------------------
n = 4;
A_true = [-0.5   0    -0.3  -0.1;
           0.4  -0.5   0.2   0  ;
           0     0.2  -0.5  -0.1;
           0.1   0.3   0    -0.5];

%--------------------------------------------------------------------------
% Simulation settings
%--------------------------------------------------------------------------
dt = 0.5;      % s
T  = 2000;     % time points
t  = (0:T-1)' * dt;

% Calcium dynamics parameters (Table 1 reference values)
tau_Ca   = 0.8;
kappa_Ca = 0.26;
R        = 0.2;
V0       = 0;
V_max    = 50;
beta_Ca  = 0.02;

% Neuronal process noise
proc_std = 0.05;

%--------------------------------------------------------------------------
% Build generative parameter struct P matching spm_fx_calcium.m's
% convention: off-diagonal P.A entries are used directly as connection
% strengths; the diagonal self-inhibition (-0.5 by construction) is
% added automatically inside spm_fx_calcium.m, so P.A's diagonal is left
% at 0 here and does not need to reproduce A_true's diagonal explicitly.
%--------------------------------------------------------------------------
P.A = A_true - diag(diag(A_true));
P.B = zeros(n,n,0);
P.C = zeros(n,0);
P.D = zeros(n,n,0);
P.tau_Ca  = log(tau_Ca/0.8);
P.kappa_Ca = log(kappa_Ca/0.26);
P.beta_Ca = zeros(n,1);
P.R       = log(R/0.2);
P.V0      = 0;

M = struct();

%--------------------------------------------------------------------------
% Simulate neuronal + calcium states (Euler integration of spm_fx_calcium)
%--------------------------------------------------------------------------
x = zeros(n,2);
x(:,2) = 0.8*(0.26*0.5 + 0.02);   % steady-state calcium at sigma(x1=0)=0.5

X = zeros(T,n,2);
u = [];
for k = 1:T
    f = spm_fx_calcium(x, u, P, M);
    f = reshape(f, n, 2);
    dx1 = f(:,1)*dt + proc_std*sqrt(dt)*randn(n,1);
    dx2 = f(:,2)*dt;
    x(:,1) = x(:,1) + dx1;
    x(:,2) = max(x(:,2) + dx2, 1e-6);
    X(k,:,:) = x;
end

%--------------------------------------------------------------------------
% Hill fluorescence observation (Table 1 reference values)
%--------------------------------------------------------------------------
Kd          = 0.2;
hill_n      = 1.0;
H_alpha_hill = 0.512;

Ca_pos  = max(squeeze(X(:,:,2)), 1e-6);   % T x n
y_clean = H_alpha_hill * (Ca_pos.^hill_n) ./ (Ca_pos.^hill_n + Kd^hill_n);

snr_db      = 15;
sig_power   = var(y_clean(:));
noise_power = sig_power / (10^(snr_db/10));
y_hill      = y_clean + sqrt(noise_power)*randn(size(y_clean));

fprintf('Simulated %d time points, %d regions. SNR = %g dB.\n', T, n, snr_db);

%--------------------------------------------------------------------------
% Assemble DCM structure. Connectivity prior is sparsity-constrained to
% the true nonzero off-diagonal entries (matching the paper's Fig. 2
% simulation design), rather than a fully-connected model.
%--------------------------------------------------------------------------
DCM_base = struct();
DCM_base.Y.dt   = dt;
DCM_base.Y.y    = y_hill;            % T x n (time bins x nodes)
DCM_base.Y.name = arrayfun(@(i) sprintf('region%d', i), 1:n, 'UniformOutput', false);

DCM_base.a = (A_true ~= 0) & ~eye(n);
DCM_base.b = zeros(n,n,0);
DCM_base.c = zeros(n,0);
DCM_base.d = zeros(n,n,0);

DCM_base.U.u    = zeros(T,1);
DCM_base.U.name = {'null'};

DCM_base.options.analysis   = 'CSD';
DCM_base.options.model      = 'Calcium';
DCM_base.options.two_state  = 0;
DCM_base.options.stochastic = 0;
DCM_base.options.induced    = 1;
DCM_base.options.maxit      = 64;
DCM_base.options.nograph    = 1;

mask = logical(DCM_base.a);

%--------------------------------------------------------------------------
% Invert under the Hill (default) observation model
%--------------------------------------------------------------------------
fprintf('\n=== Inverting with Hill observation model (default) ===\n');
DCM_hill = DCM_base;
DCM_hill = spm_dcm_calcium_csd(DCM_hill);

%--------------------------------------------------------------------------
% Invert under the linear observation model (explicit toggle)
%--------------------------------------------------------------------------
fprintf('\n=== Inverting with linear observation model ===\n');
DCM_lin = DCM_base;
DCM_lin.options.custom_g = @spm_gx_calcium_linear;
DCM_lin = spm_dcm_calcium_csd(DCM_lin);

%--------------------------------------------------------------------------
% Compare recovered connectivity to ground truth
%--------------------------------------------------------------------------
A_hill = full(DCM_hill.Ep.A);
A_lin  = full(DCM_lin.Ep.A);

r_hill    = corr(A_true(mask), A_hill(mask));
rmse_hill = sqrt(mean((A_true(mask) - A_hill(mask)).^2));

r_lin    = corr(A_true(mask), A_lin(mask));
rmse_lin = sqrt(mean((A_true(mask) - A_lin(mask)).^2));

fprintf('\n--- Hill-model estimate (r=%.3f, RMSE=%.3f) ---\n', r_hill, rmse_hill);
fprintf('--- Linear-model estimate (r=%.3f, RMSE=%.3f) ---\n', r_lin, rmse_lin);
fprintf('\nFree energy: Hill F = %.2f, Linear F = %.2f\n', DCM_hill.F, DCM_lin.F);

fprintf('\n>>>>> DEMO SCRIPT COMPLETED SUCCESSFULLY <<<<<\n');
