%% run_hidden_driver_sim.m
%
% Hidden common driver simulation for sDCM-CaI robustness check.
%
% Robustness check: adds a latent node (not included in inversion) to
% test whether spurious directed edges appear in the recovered network,
% and quantifies the resulting error.
%
% Design:
%   Ground truth: 4-node observed network (same as matched-model simulation)
%   + 1 unobserved latent driver node with coupling to k=3 of the 4 nodes.
%   Driver strength w ∈ {0, 0.05, 0.1, 0.2, 0.3, 0.5}.
%   At w=0 → identical to matched-model recovery (baseline).
%
%   For each w: N_mc = 30 Monte Carlo draws.
%   Invert sDCM-CaI on observed 4 nodes only.
%   Metrics:
%     - r(A_est, A_true)  off-diagonal correlation
%     - False positive rate (FPR): Pp>0.95 on zero connections
%     - False negative rate (FNR): Pp<0.95 on true connections
%     - Mean signed error on spurious connections
%
% OUTPUT: zebra/hidden_driver_sim_results.mat
%         figures/FigS_hidden_driver.png
%

clear; clc;

%% Config
out_mat   = './zebra/hidden_driver_sim_results.mat';
out_fig   = '../../../Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/FigS_hidden_driver.png';
w_levels  = [0, 0.05, 0.1, 0.2, 0.3, 0.5];
N_mc      = 30;
n_obs     = 4;           % observed nodes
n_total   = 5;           % 4 observed + 1 latent driver
TR        = 0.5;
T         = 2400;        % 20 min @ 2 Hz

% Ground-truth observed connectivity (same as matched-model sim)
A_obs_true = [
   -0.5   0     -0.3   -0.1
    0.4  -0.5    0.2    0
    0     0.2   -0.5   -0.1
    0.1   0.3    0     -0.5
];

% Latent driver coupling: drives nodes 1, 2, 3 (node 4 undriven)
driver_targets = [1, 2, 3];

fprintf('=== HIDDEN DRIVER SIMULATION ===\n');
fprintf('N_mc=%d, T=%d, w_levels=%s\n\n', N_mc, T, mat2str(w_levels));

%% No parallel pool — parpool unreliable on this server; parfor degrades to for

%% Pre-build reference DCM for inversion
ref_DCM            = struct();
ref_DCM.a          = ones(n_obs, n_obs);
ref_DCM.b          = zeros(n_obs, n_obs, 0);
ref_DCM.c          = zeros(n_obs, 0);
ref_DCM.d          = zeros(n_obs, n_obs, 0);
ref_DCM.Y.dt       = TR;
ref_DCM.options.analysis   = 'CSD';
ref_DCM.options.induced    = 1;
ref_DCM.options.nonlinear  = 0;
ref_DCM.options.stochastic = 0;
ref_DCM.options.maxit      = 64;
ref_DCM.options.Nmax       = 32;
ref_DCM.options.maxnodes   = 16;
ref_DCM.U.u   = zeros(n_obs, 1);
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

%% Main simulation loop
n_w   = numel(w_levels);
offdiag_mask = ~eye(n_obs);
true_mask    = A_obs_true ~= 0 & offdiag_mask;
null_mask    = A_obs_true == 0 & offdiag_mask;

r_offdiag  = nan(n_w, N_mc);
fpr        = nan(n_w, N_mc);   % false positive rate (spurious Pp>0.95)
fnr        = nan(n_w, N_mc);   % false negative rate (missed Pp<0.95)
bias_mean  = nan(n_w, N_mc);   % mean bias on true connections

for wi = 1:n_w
    w = w_levels(wi);
    fprintf('\n--- w = %.2f ---\n', w);

    % Build full 5-node A matrix (driver = node 5)
    A_full = zeros(n_total, n_total);
    A_full(1:n_obs, 1:n_obs) = A_obs_true;
    A_full(5, 5) = -0.5;              % driver self-inhibition
    for k = driver_targets
        A_full(k, 5) = w;             % driver → observed node k
    end

    Aest_all = zeros(n_obs, n_obs, N_mc);
    Pp_all   = zeros(n_obs, n_obs, N_mc);

    for mc = 1:N_mc
        rng(mc + wi*1000, 'twister');

        % Simulate 5-node system (Euler integration)
        X5 = zeros(T, n_total);
        Ca5 = zeros(T, n_total);
        tauCa = 1.6; kappa = 5.5; Qn = 0.02;
        for k = 2:T
            dX  = A_full * X5(k-1,:)' + sqrt(Qn)*randn(n_total,1);
            dCa = kappa * X5(k-1,:)' - Ca5(k-1,:)' / tauCa;
            X5(k,:)  = X5(k-1,:)  + TR * dX';
            Ca5(k,:) = Ca5(k-1,:) + TR * dCa';
        end

        % Observe only first n_obs nodes (latent driver removed)
        Y = Ca5(:, 1:n_obs);
        Y = detrend(Y);
        Y = (Y - mean(Y)) ./ (std(Y) + eps);

        % Invert sDCM-CaI on observed nodes
        DCM_mc = ref_DCM;
        DCM_mc.Y.y = Y;
        DCM_mc.M.pE = pE_ref;
        DCM_mc.M.pC = pC_ref;

        try
            DCM_est = spm_dcm_calcium_csd(DCM_mc);
            Aest_all(:,:,mc) = DCM_est.Ep.A;
            Pp_all(:,:,mc)   = DCM_est.Pp.A;
        catch
            Aest_all(:,:,mc) = NaN;
            Pp_all(:,:,mc)   = NaN;
        end

        fprintf('  w=%.2f mc=%d: F=%.1f\n', w, mc, DCM_est.F);
    end

    % Metrics per MC draw
    for mc = 1:N_mc
        A_est = Aest_all(:,:,mc);
        Pp    = Pp_all(:,:,mc);
        if all(isnan(A_est(:))), continue; end

        r_offdiag(wi,mc) = corr(A_obs_true(offdiag_mask), A_est(offdiag_mask));
        fpr(wi,mc) = mean(Pp(null_mask) > 0.95);
        fnr(wi,mc) = mean(Pp(true_mask) < 0.95);
        bias_mean(wi,mc) = mean(A_est(true_mask) - A_obs_true(true_mask));
    end

    fprintf('  w=%.2f: r=%.3f±%.3f  FPR=%.2f±%.2f  FNR=%.2f±%.2f\n', ...
        w, nanmean(r_offdiag(wi,:)), nanstd(r_offdiag(wi,:)), ...
        nanmean(fpr(wi,:)), nanstd(fpr(wi,:)), ...
        nanmean(fnr(wi,:)), nanstd(fnr(wi,:)));
end

%% Summary table
fprintf('\n%-6s  %10s  %10s  %10s  %10s\n', 'w', 'r_offdiag', 'FPR', 'FNR', 'Bias');
for wi = 1:n_w
    fprintf('%-6.2f  %10.3f  %10.3f  %10.3f  %10.3f\n', w_levels(wi), ...
        nanmean(r_offdiag(wi,:)), nanmean(fpr(wi,:)), ...
        nanmean(fnr(wi,:)), nanmean(bias_mean(wi,:)));
end

%% Figure
fig = figure('Position', [100 100 1400 460], 'Color', 'w', 'Visible', 'off');
cols_w = cool(n_w);

ax1 = subplot(1,3,1);
boxplot(r_offdiag', 'Labels', arrayfun(@(x) sprintf('%.2f',x), w_levels, 'UniformOutput', false));
hold on;
yline(nanmean(r_offdiag(1,:)), 'k--', 'LineWidth', 1.5);
ylabel('r (off-diagonal A)');
xlabel('Hidden driver strength w');
title({'Connectivity recovery'; 'vs driver strength'}, 'FontSize', 11);
grid on; box off;

ax2 = subplot(1,3,2);
r_mean = nanmean(r_offdiag, 2);
r_sd   = nanstd(r_offdiag, 0, 2);
errorbar(w_levels, r_mean, r_sd, 'ko-', 'LineWidth', 1.8, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.8]);
hold on;
fpr_mean = nanmean(fpr, 2);
plot(w_levels, fpr_mean, 'rs--', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.8 0.2 0.2]);
xlabel('Hidden driver strength w');
legend({'r_{offdiag} (recovery)', 'FPR (spurious Pp>0.95)'}, 'Location', 'best', 'FontSize', 9);
title({'Mean recovery + FPR'; 'vs driver strength'}, 'FontSize', 11);
ylim([0 1]); grid on; box off;

ax3 = subplot(1,3,3);
hold on;
for wi = 1:n_w
    plot(repmat(w_levels(wi), 1, N_mc), r_offdiag(wi,:), '.', ...
        'Color', cols_w(wi,:), 'MarkerSize', 10);
end
errorbar(w_levels, r_mean, r_sd, 'k-', 'LineWidth', 2);
xline(0.1, 'k--', 'w=0.1 threshold', 'LineWidth', 1, 'FontSize', 8);
xlabel('Hidden driver strength w');
ylabel('r (off-diagonal A)');
title({'Individual MC draws'; '(color = w level)'}, 'FontSize', 11);
grid on; box off;

% Nature-style panel labels: bold letter above top-left of each axis
for pair = {ax1,'A'; ax2,'B'; ax3,'C'}'
    ax_i = pair{1}; lbl = pair{2};
    text(ax_i, -0.08, 1.10, lbl, 'Units', 'normalized', ...
         'FontWeight', 'bold', 'FontSize', 11, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'Clipping', 'off');
end

out_dir = fileparts(out_fig);
if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
exportgraphics(fig, out_fig, 'Resolution', 600);
fprintf('\nFigure saved: %s\n', out_fig);
close(fig);

save(out_mat, 'w_levels', 'N_mc', 'r_offdiag', 'fpr', 'fnr', 'bias_mean', ...
    'A_obs_true', 'driver_targets', 'n_obs');
fprintf('Saved: %s\n', out_mat);
fprintf('\n=== Hidden driver simulation complete ===\n');
