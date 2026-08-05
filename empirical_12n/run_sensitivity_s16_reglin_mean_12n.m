%% run_sensitivity_s16_reglin_mean_12n.m
%
% Prior sensitivity analysis for RegLin-Mean 12n, Subject 16.
%
% IMPORTANT: All DCM settings are cloned from run_sc_flat_reglin_mean12n.m
% to ensure identical inversion conditions. Specifically:
%   - loads meantrace_12n.mat (same signal as reference)
%   - uses linear_obs=0 + custom_g=@spm_gx_calcium_linear (same as reference)
%   - pE.A kept at default +1/128 off-diagonal (NOT zeroed)
%   - options.pE/pC set same way as reference
%
% Two perturbations:
%   (1) Lambda scaling:  pC_pert.A = lambda * pC_ref.A  (only pC.A varies)
%   (2) Mean shift:      pE_pert.A = pE_ref.A + alpha * R .* sqrt(pC_ref.A)
%
% For lambda=1.0: identical to reference inversion -> r should be ~1.0.
%
% Output: zebra/sensitivity_12n_s16_results.mat
%         figures/Fig5.png

maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
out_mat  = fullfile(data_dir,'sensitivity_12n_s16_results.mat');
out_fig  = fullfile(home_dir,'Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/Fig5.png');

n             = 12;
N_mc          = 10;
lambda_levels = [0.25, 0.5, 1.0, 2.0, 5.0];
alpha_levels  = [0.05, 0.10, 0.25, 0.50];

%% === Clone setup from run_sc_flat_reglin_mean12n.m ===
TR = 0.5;

% Load original signal (same file as reference inversion)
d = load(fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', 16)), 'Y_mean');
signals = d.Y_mean;                          % T x 12
if size(signals,2) ~= n; signals = signals'; end
Y = detrend(signals);                        % detrend only (no demeaning/z-score)

A = ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
options.nonlinear=0; options.two_state=0; options.stochastic=0;
options.induced=1; options.centre=1; options.modality='Ca';
options.Fdcm=[0.01 0.2]; options.Nmax=32; options.maxit=256;
options.precision=log(16); options.verbose=0; options.linear_obs=0;
options.custom_g = @spm_gx_calcium_linear;

[pE_ref,pC_ref,x0] = spm_dcm_calcium_priors(A,B,C,D,options);
% Keep default +1/128 off-diagonal A prior (do NOT zero pE_ref.A)

pE_ref.alpha  = zeros(n,1);
pC_ref.alpha  = ones(n,1) * 1/16;
pE_ref.beta_y = zeros(n,1);
pC_ref.beta_y = ones(n,1) * 1/64;
pE_ref.Kd = 0; pC_ref.Kd = 1/128;
pE_ref.n  = 0; pC_ref.n  = 1/128;

pC_A = full(pC_ref.A);   % 12x12 (or 12x12x1) prior variance for A
offdiag = ~eye(n);

%% Load reference DCM to compare posteriors
ref = load(fullfile(data_dir,'subject_16_DCM_sc_flat_reglin_mean12n.mat'),'DCM_est');
DCM_ref = ref.DCM_est;
Aref    = full(DCM_ref.Ep.A);
F_ref   = DCM_ref.F;

fprintf('Reference: S16 RegLin-Mean  F=%.4f  TR=%.2f  T=%d\n', ...
    F_ref, DCM_ref.Y.dt, size(DCM_ref.Y.y,1));
fprintf('pE.A range: [%.4f %.4f] (should be ~0 to 1/128=%.4f for off-diag)\n', ...
    min(pE_ref.A(:)), max(pE_ref.A(offdiag)), 1/128);

%% Build template DCM (identical to reference inversion setup)
DCM_tpl         = struct();
DCM_tpl.a       = A;
DCM_tpl.b       = B;
DCM_tpl.c       = C;
DCM_tpl.d       = D;
DCM_tpl.Y.y     = Y;   % original detrended signal (NOT DCM_ref.Y.y)
DCM_tpl.Y.dt    = TR;
DCM_tpl.Y.Q     = [];
DCM_tpl.U.u     = [];
DCM_tpl.U.dt    = TR;
DCM_tpl.M.nograph = 1;
DCM_tpl.M.x    = x0;
DCM_tpl.options = options;
DCM_tpl.options.maxnodes = max(16, n);

%% (1) Lambda scaling — vary pC.A only, keep pE.A = default (1/128)
fprintf('\n=== Lambda scaling (N_lam=%d) ===\n', numel(lambda_levels));
A_cov    = nan(n, n, numel(lambda_levels));
F_cov    = nan(1, numel(lambda_levels));
corr_cov = nan(1, numel(lambda_levels));
sign_cov = nan(1, numel(lambda_levels));

for li = 1:numel(lambda_levels)
    lam       = lambda_levels(li);
    pC_pert   = pC_ref;
    pC_pert.A = pC_A * lam;          % scale A prior variance by lambda
    DCM_run   = DCM_tpl;
    % Only set options.pC (vary variance); options.pE uses default (pE.A=1/128)
    DCM_run.options.pE = pE_ref;     % same pE as reference (includes pE.A=1/128)
    DCM_run.options.pC = pC_pert;
    fprintf('  lambda=%.2f ...', lam);
    t0 = tic; rng(li + 100);
    try
        est  = spm_dcm_calcium_csd(DCM_run);
        Aest = full(est.Ep.A);
        Fest = est.F;
    catch ME
        fprintf(' FAILED: %s\n', ME.message);
        Aest = nan(n,n); Fest = nan;
    end
    A_cov(:,:,li) = Aest;
    F_cov(li)     = Fest;
    if ~any(isnan(Aest(:)))
        corr_cov(li) = corr(Aref(offdiag), Aest(offdiag));
        sign_cov(li) = mean(sign(Aref(offdiag)) == sign(Aest(offdiag)));
    end
    fprintf(' r=%.4f sign=%.4f dF=%.2f (%.1f min)\n', ...
        corr_cov(li), sign_cov(li), Fest - F_ref, toc(t0)/60);
end

%% (2) Alpha mean perturbation — perturb pE.A around reference prior mean (1/128)
fprintf('\n=== Alpha mean perturbation (N_mc=%d per level) ===\n', N_mc);
A_mean_mc    = nan(n, n, numel(alpha_levels), N_mc);
F_mean_all   = nan(numel(alpha_levels), N_mc);
corr_mean_all= nan(numel(alpha_levels), N_mc);
sign_mean_all= nan(numel(alpha_levels), N_mc);

for ai = 1:numel(alpha_levels)
    alp = alpha_levels(ai);
    fprintf('  alpha=%.2f :', alp);
    for mc = 1:N_mc
        rng(ai * 1000 + mc);
        R        = sign(randn(n,n));
        R(~~eye(n)) = 0;
        pE_pert   = pE_ref;
        % Perturb A around reference prior mean (pE_ref.A, off-diagonal ~1/128)
        pE_pert.A = pE_ref.A + alp * R .* sqrt(pC_A);
        pC_same   = pC_ref;   % keep pC same as reference
        DCM_run   = DCM_tpl;
        DCM_run.options.pE = pE_pert;
        DCM_run.options.pC = pC_same;
        t0 = tic;
        try
            est  = spm_dcm_calcium_csd(DCM_run);
            Aest = full(est.Ep.A);
            Fest = est.F;
        catch ME
            fprintf(' FAILED:%s', ME.message);
            Aest = nan(n,n); Fest = nan;
        end
        A_mean_mc(:,:,ai,mc) = Aest;
        F_mean_all(ai,mc)    = Fest;
        if ~any(isnan(Aest(:)))
            corr_mean_all(ai,mc) = corr(Aref(offdiag), Aest(offdiag));
            sign_mean_all(ai,mc) = mean(sign(Aref(offdiag)) == sign(Aest(offdiag)));
        end
        fprintf(' %.3f(%.0fm)', corr_mean_all(ai,mc), toc(t0)/60);
    end
    fprintf('\n  -> r=%.3f+/-%.3f  sign=%.1f+/-%.1f%%\n', ...
        nanmean(corr_mean_all(ai,:)), nanstd(corr_mean_all(ai,:)), ...
        nanmean(sign_mean_all(ai,:))*100, nanstd(sign_mean_all(ai,:))*100);
end

%% Summary
fprintf('\n=== Summary S16 RegLin-Mean 12n ===\n');
fprintf('Lambda r:    '); fprintf('%.4f ', corr_cov);    fprintf('\n');
fprintf('Lambda sign: '); fprintf('%.4f ', sign_cov);    fprintf('\n');
fprintf('Lambda range: r=%.4f--%.4f (excl. ref)\n', ...
    min(corr_cov([1 2 4 5])), max(corr_cov([1 2 4 5])));
for ai = 1:numel(alpha_levels)
    fprintf('alpha=%.2f: r=%.3f+/-%.3f  sign=%.1f+/-%.1f%%\n', ...
        alpha_levels(ai), ...
        nanmean(corr_mean_all(ai,:)), nanstd(corr_mean_all(ai,:)), ...
        nanmean(sign_mean_all(ai,:))*100, nanstd(sign_mean_all(ai,:))*100);
end

%% Save
save(out_mat, 'Aref','F_ref','A_cov','F_cov','corr_cov','sign_cov', ...
    'A_mean_mc','F_mean_all','corr_mean_all','sign_mean_all', ...
    'lambda_levels','alpha_levels','N_mc','n','-v7.3');
fprintf('Saved: %s\n', out_mat);

%% Figure (4 panels)
fig = figure('Position',[100 100 3200 800],'Color','w','Visible','off');

ax1 = subplot(1,4,1);
plot(lambda_levels, corr_cov, 'ko-','MarkerFaceColor','k','LineWidth',2,'MarkerSize',8);
xlabel('\lambda','FontSize',12); ylabel('r with reference','FontSize',12);
title('(A) Prior variance scaling','FontSize',12);
ylim([0 1.05]); grid on; box off; xline(1,'k--','LineWidth',1);

ax2 = subplot(1,4,2);
plot(lambda_levels, sign_cov*100, 'bs-','MarkerFaceColor','b','LineWidth',2,'MarkerSize',8);
xlabel('\lambda','FontSize',12); ylabel('Sign consistency (%)','FontSize',12);
title('(B) Sign (\lambda)','FontSize',12);
ylim([0 105]); grid on; box off; xline(1,'k--','LineWidth',1);

ax3 = subplot(1,4,3);
r_mean = nanmean(corr_mean_all, 2);
r_std  = nanstd(corr_mean_all, 0, 2);
errorbar(alpha_levels, r_mean, r_std, 'ro-','MarkerFaceColor','r','LineWidth',2,'MarkerSize',8,'CapSize',6);
xlabel('\alpha','FontSize',12); ylabel('r with reference','FontSize',12);
title('(C) Prior mean perturbation','FontSize',12);
ylim([0 1.05]); grid on; box off;

ax4 = subplot(1,4,4);
s_mean = nanmean(sign_mean_all, 2)*100;
s_std  = nanstd(sign_mean_all, 0, 2)*100;
errorbar(alpha_levels, s_mean, s_std, 'ms-','MarkerFaceColor','m','LineWidth',2,'MarkerSize',8,'CapSize',6);
xlabel('\alpha','FontSize',12); ylabel('Sign consistency (%)','FontSize',12);
title('(D) Sign (\alpha)','FontSize',12);
ylim([0 105]); grid on; box off;

sgtitle(sprintf('S16 RegLin-Mean 12n prior sensitivity (N_{mc}=%d, T=%d)', ...
    N_mc, size(Y,1)), 'FontSize',13,'FontWeight','bold');

exportgraphics(fig, out_fig, 'Resolution',300);
fprintf('Figure saved: %s\n', out_fig);
close(fig);
fprintf('Done.\n');
