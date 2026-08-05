%% run_mc_reglin_mean_12n.m
%
% Monte Carlo parameter recovery for RegLin-Mean 12n.
% Calibrated from real S12 RegLin-Mean DCM (options, TR, T).
% Same SDE simulation as run_mc_simulation_12node.m, but explicitly
% tied to the RegLin-Mean model selected by BMS.
%
% OUTPUT: zebra/mc_reglin_mean_12n_results.mat
%         figures/FigS_mc_reglin_mean_12n.png

clear; clc;
maxNumCompThreads(1);
rng(42);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
out_mat  = fullfile(data_dir,'mc_reglin_mean_12n_results.mat');
out_fig  = fullfile(home_dir,'Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/FigS_mc_reglin_mean_12n.png');
N_mc     = 20;
n        = 12;

node_names = {'lTeO','lTh','lP','lPT','lHb','lpRF', ...
              'rTeO','rTh','rP','rPT','rHb','rpRF'};

%% Load RegLin-Mean reference DCM (S12) for options/TR/T
ref_file = fullfile(data_dir, 'subject_12_DCM_sc_flat_reglin_mean12n.mat');
ref      = load(ref_file, 'DCM_est');
ref_DCM0 = ref.DCM_est;
TR       = ref_DCM0.Y.dt;
T        = size(ref_DCM0.Y.y, 1);
options  = ref_DCM0.options;
options.verbose  = 0;
options.maxit    = 64;
options.Nmax     = 32;
options.maxnodes = max(16, n);
fprintf('RegLin-Mean reference: S12  T=%d  TR=%.2f  linear_obs=%d\n', T, TR, options.linear_obs);

%% Build sparse stable A_true (~20% off-diagonal density)
rng(0);
A_true = -0.5 * eye(n);
n_conn = 26;
offdiag_ij = zeros(0,2);
while size(offdiag_ij,1) < n_conn
    i = randi(n); j = randi(n);
    if i ~= j && ~any(offdiag_ij(:,1)==i & offdiag_ij(:,2)==j)
        offdiag_ij(end+1,:) = [i j];
    end
end
for k = 1:n_conn
    i = offdiag_ij(k,1); j = offdiag_ij(k,2);
    A_true(i,j) = 0.3 * randn();
end
for attempt = 1:20
    ev = real(eig(A_true));
    if max(ev) < -0.05, break; end
    A_true = A_true * 0.9;
end
fprintf('A_true: max Re(eig)=%.4f  non-zero off-diag=%d\n\n', max(real(eig(A_true))), nnz(A_true & ~eye(n)));

%% Build reference DCM structure from RegLin-Mean options
A_tpl = 0.001 * ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
[pE_ref, pC_ref, x0_ref] = spm_dcm_calcium_priors(A_tpl,B,C,D,options);
pE_ref.A = pE_ref.A * 0;

ref_DCM            = struct();
ref_DCM.a          = double(A_true ~= 0);
ref_DCM.b          = B; ref_DCM.c = C; ref_DCM.d = D;
ref_DCM.Y.dt       = TR;
ref_DCM.Y.Q        = [];
ref_DCM.U.u        = []; ref_DCM.U.dt = TR;
ref_DCM.M.pE       = pE_ref;
ref_DCM.M.pC       = pC_ref;
ref_DCM.M.x        = x0_ref;
ref_DCM.M.n        = size(x0_ref,1);
ref_DCM.M.m        = 0;
ref_DCM.M.l        = size(x0_ref,2);
ref_DCM.M.nograph  = 1;
ref_DCM.options    = options;

pE_true   = pE_ref;
pE_true.A = A_true;

%% Monte Carlo recovery
Aest_all = nan(n,n,N_mc);
Pp_all   = nan(n,n,N_mc);
F_all    = nan(1,N_mc);

for mc = 1:N_mc
    rng(mc, 'twister');
    fprintf('--- MC rep %d/%d ---\n', mc, N_mc);

    % Simulate calcium dynamics from A_true
    X  = zeros(T,n); Ca = zeros(T,n);
    tauCa = 1.6; kappa = 5.5; Qn = 0.02;
    for k = 2:T
        dX  = A_true * X(k-1,:)' + sqrt(Qn)*randn(n,1);
        dCa = kappa * X(k-1,:)' - Ca(k-1,:)' / tauCa;
        X(k,:)  = X(k-1,:)  + TR * dX';
        Ca(k,:) = Ca(k-1,:) + TR * dCa';
    end
    Y = detrend(Ca);
    Y = (Y - mean(Y)) ./ (std(Y) + eps);

    DCM_mc         = ref_DCM;
    DCM_mc.Y.y     = Y;
    DCM_mc.M.pE    = pE_true;
    DCM_mc.name    = sprintf('MC_reglin_mean_12n_rep%02d', mc);

    t0 = tic;
    try
        DCM_est = spm_dcm_calcium_csd(DCM_mc);
        Aest_all(:,:,mc) = DCM_est.Ep.A;
        Pp_all(:,:,mc)   = DCM_est.Pp.A;
        F_all(mc)        = DCM_est.F;
        fprintf('  F=%.2f  (%.1f min)\n', DCM_est.F, toc(t0)/60);
    catch ME
        fprintf('  FAILED: %s\n', ME.message);
    end
end

%% Metrics
offdiag_mask = ~eye(n);
true_mask    = A_true ~= 0 & offdiag_mask;

r_vals    = nan(1,N_mc);
rmse_vals = nan(1,N_mc);
sign_acc  = nan(1,N_mc);
auc_vals  = nan(1,N_mc);

for mc = 1:N_mc
    if all(isnan(Aest_all(:,:,mc))), continue; end
    A_est = Aest_all(:,:,mc);
    r_vals(mc)    = corr(A_true(offdiag_mask), A_est(offdiag_mask));
    rmse_vals(mc) = sqrt(mean((A_true(true_mask) - A_est(true_mask)).^2));
    sign_acc(mc)  = mean(sign(A_true(offdiag_mask)) == sign(A_est(offdiag_mask)));
    try
        [~,~,~,auc_vals(mc)] = perfcurve(true_mask(offdiag_mask), Pp_all(offdiag_mask + n^2*(mc-1)), 1);
    catch, end
end

fprintf('\n=== MC RegLin-Mean 12n (N=%d, T=%d) ===\n', N_mc, T);
fprintf('r (off-diagonal):  %.3f +/- %.3f\n', nanmean(r_vals), nanstd(r_vals));
fprintf('RMSE (true edges): %.4f +/- %.4f\n', nanmean(rmse_vals), nanstd(rmse_vals));
fprintf('Sign accuracy:     %.1f +/- %.1f%%\n', nanmean(sign_acc)*100, nanstd(sign_acc)*100);
fprintf('AUC:               %.3f +/- %.3f\n', nanmean(auc_vals), nanstd(auc_vals));
fprintf('Mean F:            %.2f +/- %.2f\n', nanmean(F_all), nanstd(F_all));

%% Figure
Amean = nanmean(Aest_all,3);
fig = figure('Position',[100 100 2400 800],'Color','w','Visible','off');

ax1 = subplot(1,3,1);
scatter(A_true(offdiag_mask), Amean(offdiag_mask), 80, [0.2 0.5 0.8], 'filled', 'MarkerFaceAlpha',0.7);
hold on;
lims = [min([A_true(offdiag_mask); Amean(offdiag_mask)]) max([A_true(offdiag_mask); Amean(offdiag_mask)])];
plot(lims, lims, 'k--', 'LineWidth',1.5);
xlabel('A_{true}','FontSize',12); ylabel('A_{est} (mean)','FontSize',12);
title(sprintf('(A) Recovery\nr=%.3f±%.3f', nanmean(r_vals), nanstd(r_vals)),'FontSize',12);
axis square; grid on; box off;

ax2 = subplot(1,3,2);
bar(1:N_mc, r_vals, 0.7, 'FaceColor',[0.2 0.5 0.8],'EdgeColor','none');
yline(nanmean(r_vals),'r--','LineWidth',2);
xlabel('MC run','FontSize',12); ylabel('r','FontSize',12);
title(sprintf('(B) Per-run r\n(mean=%.3f)', nanmean(r_vals)),'FontSize',12);
ylim([-0.1 1]); grid on; box off;

ax3 = subplot(1,3,3);
mx = max(abs(A_true(:)));
imagesc(A_true,[-mx mx]); colormap(ax3, redblue(256));
colorbar; axis square;
set(ax3,'XTick',1:n,'XTickLabel',node_names,'XTickLabelRotation',45,...
    'YTick',1:n,'YTickLabel',node_names,'FontSize',8);
title('(C) A_{true}','FontSize',12);

sgtitle(sprintf('RegLin-Mean 12n face validity (N_{mc}=%d, T=%d, TR=%.1fs)', N_mc, T, TR),...
    'FontSize',13,'FontWeight','bold');

exportgraphics(fig, out_fig, 'Resolution',300);
fprintf('Figure saved: %s\n', out_fig);
close(fig);

save(out_mat,'A_true','Aest_all','Pp_all','F_all',...
    'r_vals','rmse_vals','sign_acc','auc_vals','node_names','N_mc','T','TR');
fprintf('Saved: %s\nDone.\n', out_mat);

function cmap = redblue(n)
half = floor(n/2);
r = [linspace(0.1,1,half), linspace(1,0.7,n-half)];
g = [linspace(0.3,1,half), linspace(1,0.1,n-half)];
b = [linspace(0.8,1,half), linspace(1,0.1,n-half)];
cmap = [r(:),g(:),b(:)];
end
