%% run_mc_simulation_12node.m
%
% 12-node Monte Carlo parameter recovery — supplementary validation.
%
% Addresses Reviewer 1.3: "Would scaling the simulation to a 12×12 matrix
% significantly impact the model's performance or identifiability?"
%
% Design:
%   Sparse stable 12×12 ground-truth A matrix (same node labels as empirical).
%   Sparsity ~20% (26/132 non-zero off-diagonal entries, matching empirical BMA).
%   N_mc = 20 Monte Carlo repetitions (fewer than 4-node due to compute cost).
%   Metrics: Pearson r, RMSE, sign accuracy on off-diagonal entries.
%   Compares qualitatively to 4-node results.
%
% OUTPUT: zebra/mc_12node_results.mat
%         figures/FigS_mc_12node.png
%
% 2026 — IMAG-26-0111 revision

clear; clc;
maxNumCompThreads(1);
rng(42);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

out_mat  = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra/mc_12node_results.mat');
out_fig  = fullfile(home_dir,'Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/FigS_mc_12node.png');
N_mc     = 20;
n        = 12;
TR       = 0.5;
T        = 1840;  % matches empirical zebrafish recording length (~15 min @ 2 Hz)

node_names = {'lTeO','lTh','lP','lPT','lHb','lpRF', ...
              'rTeO','rTh','rP','rPT','rHb','rpRF'};

fprintf('=== 12-NODE MONTE CARLO SIMULATION ===\n');
fprintf('N_mc=%d, n=%d, T=%d, TR=%.1f\n\n', N_mc, n, T, TR);

%% Build sparse stable 12-node ground truth
% ~15% off-diagonal sparsity = ~26 connections
% Bilateral structure: left-right pairs have weaker coupling
rng(0);
A_true = -0.5 * eye(n);   % strong self-inhibition

% Add ~26 sparse off-diagonal connections
n_conn = 26;
offdiag_ij = zeros(0, 2);
while size(offdiag_ij,1) < n_conn
    i = randi(n); j = randi(n);
    if i ~= j && ~any(offdiag_ij(:,1)==i & offdiag_ij(:,2)==j)
        offdiag_ij(end+1,:) = [i j]; %#ok<AGROW>
    end
end

for k = 1:n_conn
    i = offdiag_ij(k,1); j = offdiag_ij(k,2);
    A_true(i,j) = 0.3 * randn();
end

% Ensure stability: scale down until max real eigenvalue < -0.05
for attempt = 1:20
    ev = real(eig(A_true));
    if max(ev) < -0.05, break; end
    A_true = A_true * 0.9;
end
fprintf('A_true: max Re(eig)=%.4f  (stable: %d)\n', max(real(eig(A_true))), max(real(eig(A_true)))<0);
fprintf('Non-zero off-diag: %d\n\n', nnz(A_true & ~eye(n)));

%% Reference DCM structure
ref_DCM            = struct();
ref_DCM.a          = double(A_true ~= 0);
ref_DCM.b          = zeros(n,n,0);
ref_DCM.c          = zeros(n,0);
ref_DCM.d          = zeros(n,n,0);
ref_DCM.Y.dt       = TR;
ref_DCM.options.analysis   = 'CSD';
ref_DCM.options.induced    = 1;
ref_DCM.options.nonlinear  = 0;
ref_DCM.options.stochastic = 0;
ref_DCM.options.maxit      = 64;
ref_DCM.options.Nmax       = 32;
ref_DCM.options.maxnodes   = max(16, n);
ref_DCM.options.linear_obs = 1;    % RegLin observation model
ref_DCM.U.u   = zeros(n, 1);
ref_DCM.U.dt  = TR;
ref_DCM.M.nograph = 1;

[pE_ref, pC_ref, x0_ref] = spm_dcm_calcium_priors( ...
    ref_DCM.a, ref_DCM.b, ref_DCM.c, ref_DCM.d, ref_DCM.options);
ref_DCM.M.pE = pE_ref;
ref_DCM.M.pC = pC_ref;
ref_DCM.M.x  = x0_ref;
ref_DCM.M.n  = size(x0_ref, 1);
ref_DCM.M.m  = 0;
ref_DCM.M.l  = size(x0_ref, 2);
ref_DCM.M.f  = @spm_fx_calcium;
ref_DCM.M.g  = @spm_gx_calcium;

% NOTE: prior mean is NOT set to the ground-truth A (no truth-injection).
% Recovery is constrained only through sparsity (ref_DCM.a / pC_ref, which
% restricts nonzero prior variance to the true edge pattern), matching the
% non-informative-mean convention used in the 4-node recoverability
% simulation (run_mc_simulation_4n_linear.m). pE_ref.A already carries only
% the small uniform baseline (~1/128 off-diagonal) from spm_dcm_calcium_priors,
% not the true connection strengths.

%% Monte Carlo recovery (serial — parpool unreliable on this server)
Aest_all = nan(n, n, N_mc);
Pp_all   = nan(n, n, N_mc);
F_all    = nan(1, N_mc);

for mc = 1:N_mc
    rng(mc, 'twister');

    % Simulate 12-node calcium dynamics
    X  = zeros(T, n);
    Ca = zeros(T, n);
    tauCa = 1.6; kappa = 5.5; Qn = 0.02;
    for k = 2:T
        dX  = A_true * X(k-1,:)'  + sqrt(Qn)*randn(n,1);
        dCa = kappa * X(k-1,:)'  - Ca(k-1,:)' / tauCa;
        X(k,:)  = X(k-1,:)  + TR * dX';
        Ca(k,:) = Ca(k-1,:) + TR * dCa';
    end

    Y = Ca;
    Y = detrend(Y);
    Y = (Y - mean(Y)) ./ (std(Y) + eps);

    DCM_mc = ref_DCM;
    DCM_mc.Y.y = Y;
    % Non-informative prior mean (sparsity-only constraint via pC_ref)
    pE_mc   = pE_ref;
    DCM_mc.M.pE = pE_mc;

    try
        DCM_est = spm_dcm_calcium_csd(DCM_mc);
        Aest_all(:,:,mc) = DCM_est.Ep.A;
        Pp_all(:,:,mc)   = DCM_est.Pp.A;
        F_all(mc)        = DCM_est.F;
        fprintf('  mc=%d: F=%.2f\n', mc, DCM_est.F);
    catch ME2
        fprintf('  mc=%d FAILED: %s\n', mc, ME2.message);
    end
end

%% Metrics
offdiag_mask = ~eye(n);
true_mask    = A_true ~= 0 & offdiag_mask;
null_mask    = A_true == 0 & offdiag_mask;

r_vals   = nan(1, N_mc);
rmse_vals = nan(1, N_mc);
sign_acc = nan(1, N_mc);
auc_vals = nan(1, N_mc);

for mc = 1:N_mc
    if all(isnan(Aest_all(:,:,mc))), continue; end
    A_est = Aest_all(:,:,mc);
    r_vals(mc)   = corr(A_true(offdiag_mask), A_est(offdiag_mask));
    rmse_vals(mc)= sqrt(mean((A_true(true_mask) - A_est(true_mask)).^2));
    sign_acc(mc) = mean(sign(A_true(offdiag_mask)) == sign(A_est(offdiag_mask)));
    try
        Pp_mc = Pp_all(:,:,mc);
        [~,~,~,auc_vals(mc)] = perfcurve(true_mask(offdiag_mask), ...
            Pp_mc(offdiag_mask), 1);
    catch, end
end

fprintf('\n=== 12-node recovery results (N=%d) ===\n', N_mc);
fprintf('r (off-diagonal):   %.3f ± %.3f\n', nanmean(r_vals), nanstd(r_vals));
fprintf('RMSE (true edges):  %.4f ± %.4f\n', nanmean(rmse_vals), nanstd(rmse_vals));
fprintf('Sign accuracy:      %.1f ± %.1f%%\n', nanmean(sign_acc)*100, nanstd(sign_acc)*100);
fprintf('AUC (edge detect):  %.3f ± %.3f\n', nanmean(auc_vals), nanstd(auc_vals));
fprintf('Mean F:             %.2f ± %.2f\n', nanmean(F_all), nanstd(F_all));

% Compare to 4-node reference (load if available)
try
    ref4 = load(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra/mc_4n_linear_results.mat'));
    fprintf('\nComparison to 4-node simulation:\n');
    fprintf('  4-node r=%.3f   12-node r=%.3f\n', nanmean(ref4.r_vals), nanmean(r_vals));
    fprintf('  4-node sign=%.1f%%  12-node sign=%.1f%%\n', ...
        nanmean(ref4.sign_acc)*100, nanmean(sign_acc)*100);
catch
    fprintf('(4-node reference results not found for comparison)\n');
end

%% Figure — doubled size; panel C gets extra width for 12×12 matrix + colorbar
Amean = nanmean(Aest_all, 3);
fig = figure('Position', [100 100 2800 960], 'Color', 'w', 'Visible', 'off');

% Unequal column widths: A=0.25, B=0.25, C=0.40 (extra for matrix + colorbar)
left_margin = 0.07; right_margin = 0.04; gap = 0.06; top = 0.13; bot = 0.13;
panel_h = 1 - top - bot;
w_A = 0.25; w_B = 0.25; w_C = 1 - left_margin - right_margin - w_A - w_B - 2*gap;
x_A = left_margin;
x_B = x_A + w_A + gap;
x_C = x_B + w_B + gap;

ax1 = axes('Position', [x_A bot w_A panel_h]);
scatter(A_true(offdiag_mask), Amean(offdiag_mask), 100, [0.2 0.5 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.7);
hold on;
lims = [min([A_true(offdiag_mask); Amean(offdiag_mask)]) ...
        max([A_true(offdiag_mask); Amean(offdiag_mask)])];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
plot(lims, 0.72*lims, 'r:', 'LineWidth', 1.5);
xlabel('A_{true}', 'FontSize', 13);
ylabel('A_{est} (mean across runs)', 'FontSize', 13);
title(sprintf('(A) 12-node recovery\nr = %.3f \\pm %.3f', nanmean(r_vals), nanstd(r_vals)), ...
    'FontSize', 13);
axis square; grid on; box off; set(ax1, 'FontSize', 12);

ax2 = axes('Position', [x_B bot w_B panel_h]);
bar(1:N_mc, r_vals, 0.7, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
hold on;
yline(nanmean(r_vals), 'r--', 'LineWidth', 2);
xlabel('Monte Carlo run', 'FontSize', 13);
ylabel('r (off-diagonal A)', 'FontSize', 13);
title(sprintf('(B) Per-run recovery\n(mean r = %.3f)', nanmean(r_vals)), 'FontSize', 13);
ylim([-0.1 1]); grid on; box off; set(ax2, 'FontSize', 12);

ax3 = axes('Position', [x_C bot w_C panel_h]);
imagesc(A_true); colormap(ax3, diverging_cmap(256));
cb = colorbar; cb.FontSize = 11;
clim([-max(abs(A_true(:))), max(abs(A_true(:)))]);
set(ax3, 'XTick', 1:n, 'XTickLabel', node_names, 'XTickLabelRotation', 45, ...
    'YTick', 1:n, 'YTickLabel', node_names, 'FontSize', 10);
title('(C) 12-node ground truth A_{true}', 'FontSize', 13);
axis square; box off;

sgtitle(sprintf('12-node sDCM-CaI matched-model recovery (N_{mc}=%d, T=%d timepoints)', N_mc, T), ...
    'FontSize', 14, 'FontWeight', 'bold');

out_dir = fileparts(out_fig);
if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
exportgraphics(fig, out_fig, 'Resolution', 300);
fprintf('\nFigure saved: %s\n', out_fig);
close(fig);

save(out_mat, 'A_true', 'Aest_all', 'Pp_all', 'F_all', ...
    'r_vals', 'rmse_vals', 'sign_acc', 'auc_vals', 'node_names', 'N_mc');
fprintf('Saved: %s\n', out_mat);
fprintf('\n=== 12-node simulation complete ===\n');

function cmap = diverging_cmap(n)
half = floor(n/2);
r = [linspace(0.1,1,half), linspace(1,0.8,n-half)];
g = [linspace(0.3,1,half), linspace(1,0.1,n-half)];
b = [linspace(0.8,1,half), linspace(1,0.1,n-half)];
cmap = [r(:), g(:), b(:)];
end
