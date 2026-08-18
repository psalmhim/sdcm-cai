%% misspecified_model_sim.m
%
% Model-misspecification recovery test for spectral DCM-CaI (single subject)
%
% Design:
%   TRUE generative model : Linear SDE + direct linear observation (no calcium)
%                           dx = A_phys * x * dt_sim + sigma * dW
%                           y  = x + meas_noise      (no Hill, no Ca kinetics)
%   INVERTING model       : sDCM-CaI (assumes calcium kinetics + Hill obs)
%
% Key design choice: A_true is a synthetic sparse 12-node network with
% meaningful off-diagonal coupling (0.1–0.3 range), same sparsity as
% the 12-node MC simulation elsewhere in this repo. Individual-level S16 Ep.A is NOT used
% because VB-Laplace shrinks individual posteriors toward zero (off-diagonal
% ≈ 0.002–0.005), which would trivialize the recovery test.
%
% Stability: sub-stepped Euler (dt_sim = dt/10) avoids instability
% from large diagonal self-inhibition.
%
% Output: r, sign consistency, RMSE for off-diagonal A
%

clear; clc;
rng(42);

home_dir = getenv('HOME');
addpath(fullfile(home_dir, 'Dropbox/matlabwork/mnet0.92/dcmcai'));
addpath(genpath(fullfile(home_dir, 'Dropbox/matlabwork/spm25')));

%% 1. Load S16 DCM to get inversion options and region names
dcm_file = fullfile(home_dir, 'Dropbox/matlabwork/mnet0.92/dcmcai/DCM_cai_pc1active_subject_16.mat');
fprintf('Loading %s for inversion settings...\n', dcm_file);
D = load(dcm_file);
DCM_S16 = D.DCM;
n  = size(DCM_S16.Ep.A, 1);   % 12
dt = DCM_S16.M.dt;             % 0.5 s (2 Hz)

%% 2. Construct synthetic sparse 12-node A_true (physical connectivity)
%
% Sparse design: ~20% off-diagonal density (≈26/132 connections nonzero)
% Coupling magnitudes: excitatory 0.10–0.25, inhibitory -0.10 to -0.20
% Diagonal: -0.5 (moderate self-inhibition, stable)
%
% Same sparse pattern as the 12-node MC simulation for comparability.

A_true = diag(-0.5 * ones(1,n));   % diagonal self-inhibition

% Sparse connectivity seed (reproducible, ~20% density)
rng(7);
density   = 0.20;
n_nonzero = round(density * n*(n-1));   % ~26 off-diagonal nonzero

% Random sparse coupling
idx_offdiag = find(~eye(n));
perm        = randperm(numel(idx_offdiag));
active_idx  = idx_offdiag(perm(1:n_nonzero));
for k = 1:n_nonzero
    sign_k = sign(randn);
    mag_k  = 0.10 + 0.15*rand;    % uniform 0.10 to 0.25
    A_true(active_idx(k)) = sign_k * mag_k;
end

% Ensure stability (max real eigenvalue < 0)
ev = max(real(eig(A_true)));
fprintf('Max real eigenvalue of A_true: %.4f  ', ev);
if ev >= 0
    A_true = A_true - (ev + 0.05)*eye(n);
    fprintf('[stabilized -> %.4f]', max(real(eig(A_true))));
end
fprintf('\n');
fprintf('Off-diagonal nonzero: %d / %d\n', n_nonzero, n*(n-1));
fprintf('A_true off-diag range (nonzero): %.3f to %.3f\n', ...
    min(A_true(active_idx)), max(A_true(active_idx)));

%% 3. Generate signals from LINEAR SDE with sub-stepped Euler
%    True model: dx = A_true*x*dt_sim + sigma*dW, y = x + noise
%    Sub-step with dt_sim = dt/10 to stabilize Euler for |A_diag| ≤ 0.5

T        = 1600;       % 800 s at 2 Hz (longer than S16's 672 for better identif.)
n_sub    = 10;         % Euler sub-steps per observed sample
dt_sim   = dt/n_sub;   % 0.05 s sub-step
sigma    = 0.10;       % endogenous noise
meas_sd  = 0.02;       % measurement noise

x = zeros(n, 1);
Y = zeros(T, n);

rng(42);
for t = 1:T
    for s = 1:n_sub
        x = x + dt_sim * (A_true * x) + sqrt(dt_sim)*sigma*randn(n,1);
    end
    Y(t,:) = x' + meas_sd*randn(1,n);
end

fprintf('\nSignal statistics (linear SDE, no calcium):\n');
fprintf('  Mean abs(y): %.4f\n', mean(abs(Y(:))));
fprintf('  Std(y):      %.4f\n', std(Y(:)));

if any(isnan(Y(:))) || any(isinf(Y(:)))
    error('Signal is NaN/Inf — system unstable. Check A_true eigenvalues.');
end

%% 4. Build DCM structure for sDCM-CaI inversion
DCM_mis           = struct();
DCM_mis.Y.y       = Y;
DCM_mis.Y.dt      = dt;
% Region labels
region_labels = cell(1,n);
for i=1:n, region_labels{i} = sprintf('R%02d',i); end
DCM_mis.Y.name    = region_labels;
DCM_mis.v         = T;
DCM_mis.n         = n;

DCM_mis.a         = ones(n,n);    % fully connected prior (same as empirical)
DCM_mis.b         = zeros(n,n,0);
DCM_mis.c         = zeros(n,0);
DCM_mis.d         = zeros(n,n,0);

DCM_mis.U.u       = zeros(T,1);
DCM_mis.U.name    = {'null'};

% Copy inversion options from S16 (frequency band, order, etc.)
DCM_mis.options          = DCM_S16.options;
DCM_mis.options.nograph  = 1;

%% 5. Priors (same as main empirical analysis)
[pE, pC, x0] = spm_dcm_calcium_priors(DCM_mis.a, DCM_mis.b, ...
                                        DCM_mis.c, DCM_mis.d, DCM_mis.options);
DCM_mis.M.f  = DCM_S16.M.f;
DCM_mis.M.g  = DCM_S16.M.g;
DCM_mis.M.x  = x0;
DCM_mis.M.dt = dt;
DCM_mis.M.N  = DCM_S16.M.N;
DCM_mis.options.pE = pE;
DCM_mis.options.pC = pC;

%% 6. Invert with sDCM-CaI (CALCIUM model on LINEAR data — misspecification)
fprintf('\nInverting linear-SDE data with sDCM-CaI (misspecified)...\n');
t_start = tic;
DCM_est = spm_dcm_calcium_csd(DCM_mis);
elapsed = toc(t_start)/60;
fprintf('Inversion completed in %.1f min\n', elapsed);

%% 7. Evaluate recovery (off-diagonal only)
offdiag  = logical(~eye(n));
nonzero  = false(n,n);
nonzero(active_idx) = true;

% Compare recovered Ep.A (off-diagonal) vs A_true (off-diagonal, physical)
% Note: off-diagonal Ep.A ≈ physical coupling in DCM parameterization
v_true_all  = A_true(offdiag);             % true off-diagonal (all 132)
v_rec_all   = DCM_est.Ep.A(offdiag);      % recovered (all 132)

v_true_nz   = A_true(nonzero & offdiag);   % true nonzero only
v_rec_nz    = DCM_est.Ep.A(nonzero & offdiag);

r_all       = corr(v_true_all, v_rec_all);
r_nonzero   = corr(v_true_nz, v_rec_nz);
sign_all    = mean(sign(v_true_all) == sign(v_rec_all));
sign_nonzero = mean(sign(v_true_nz) == sign(v_rec_nz));
rmse_all    = sqrt(mean((v_true_all - v_rec_all).^2));

fprintf('\n========== MISSPECIFIED MODEL RECOVERY RESULTS ==========\n');
fprintf('TRUE model:      12-node linear SDE, direct obs (no Ca kinetics)\n');
fprintf('INVERTING model: sDCM-CaI (calcium kinetics + Hill obs model)\n');
fprintf('A_true:          sparse 12-node, %d nonzero off-diagonal (%.0f%%)\n', ...
    n_nonzero, density*100);
fprintf('----------------------------------------------------------\n');
fprintf('r (all 132 off-diagonal):        %.3f\n', r_all);
fprintf('r (nonzero %d off-diagonal only): %.3f\n', n_nonzero, r_nonzero);
fprintf('Sign cons. (all 132):            %.1f%%\n', sign_all*100);
fprintf('Sign cons. (nonzero %d only):     %.1f%%\n', n_nonzero, sign_nonzero*100);
fprintf('RMSE (all 132 off-diagonal):     %.4f\n', rmse_all);
fprintf('==========================================================\n');

%% 8. Save
save('zebra/misspecified_model_results.mat', ...
    'r_all', 'r_nonzero', 'sign_all', 'sign_nonzero', 'rmse_all', ...
    'A_true', 'DCM_est', 'n_nonzero', 'density', 'elapsed');
fprintf('Results saved to zebra/misspecified_model_results.mat\n');

%% 9. Figure
fig = figure('Visible','off');
scatter(v_true_all, v_rec_all, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(v_true_nz, v_rec_nz, 60, 'r', 'filled');
lims = [min([v_true_all;v_rec_all])-0.05, max([v_true_all;v_rec_all])+0.05];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
xlabel('True A_{ij} (linear SDE)');
ylabel('Recovered A_{ij} (sDCM-CaI, misspecified)');
title(sprintf('Misspecified model recovery  r=%.3f (all)  sign=%.0f%%', ...
    r_all, sign_all*100));
legend({'All off-diagonal (n=132)','Nonzero connections','Identity'},'Location','nw');
grid on;
saveas(fig, 'zebra/misspecified_model_scatter.png');
fprintf('Figure saved.\n');
