%% run_sensitivity_s16_merge.m
% Merge alpha_a + alpha_b results with hard-coded lambda results,
% save final sensitivity_12n_s16_results.mat and regenerate Fig5.png.
% Run after BOTH alpha_a and alpha_b are complete.

home_dir = getenv('HOME');
data_dir = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
out_mat  = fullfile(data_dir,'sensitivity_12n_s16_results.mat');
out_fig  = fullfile(home_dir,'Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/Fig5.png');

n = 12;
alpha_levels_all = [0.05, 0.10, 0.25, 0.50];
lambda_levels    = [0.25, 0.5, 1.0, 2.0, 5.0];
N_mc = 10;

%% Hard-coded lambda results (from run_sensitivity_s16_reglin_mean_12n.m, completed 2026-07-21)
corr_cov = [0.9642, 0.9916, 1.0000, 0.9935, 0.9925];
sign_cov = [0.9015, 0.9545, 1.0000, 0.9545, 0.9621];
F_ref    = -3632.7384;

%% Load alpha results
% alpha_a: alpha=[0.05,0.10], rows 1-2, fields corr_mean_all (2x10)
% alpha_b025: alpha=0.25, fields corr_mc (1x10)
% alpha_b050: alpha=0.50, fields corr_mc (1x10)
a    = load(fullfile(data_dir,'sensitivity_12n_s16_alpha_a.mat'));
b025 = load(fullfile(data_dir,'sensitivity_12n_s16_alpha_b025.mat'));
b050 = load(fullfile(data_dir,'sensitivity_12n_s16_alpha_b050.mat'));

corr_mean_all = [a.corr_mean_all; b025.corr_mc; b050.corr_mc];  % 4 x 10
sign_mean_all = [a.sign_mean_all; b025.sign_mc; b050.sign_mc];  % 4 x 10

fprintf('=== Merged sensitivity S16 RegLin-Mean 12n ===\n');
fprintf('Lambda r:    '); fprintf('%.4f ', corr_cov); fprintf('\n');
fprintf('Lambda sign: '); fprintf('%.4f ', sign_cov); fprintf('\n');
fprintf('Lambda range (excl. ref): r=%.4f--%.4f\n', min(corr_cov([1 2 4 5])), max(corr_cov([1 2 4 5])));
for ai = 1:4
    fprintf('alpha=%.2f: r=%.3f+/-%.3f  sign=%.1f+/-%.1f%%\n', ...
        alpha_levels_all(ai), ...
        nanmean(corr_mean_all(ai,:)), nanstd(corr_mean_all(ai,:)), ...
        nanmean(sign_mean_all(ai,:))*100, nanstd(sign_mean_all(ai,:))*100);
end

save(out_mat, 'F_ref','corr_cov','sign_cov','corr_mean_all','sign_mean_all', ...
    'lambda_levels','alpha_levels_all','N_mc','n','-v7.3');
fprintf('Saved: %s\n', out_mat);

%% Figure
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
fig = figure('Position',[100 100 3200 800],'Color','w','Visible','off');

subplot(1,4,1);
plot(lambda_levels, corr_cov, 'ko-','MarkerFaceColor','k','LineWidth',2,'MarkerSize',8);
xlabel('\lambda','FontSize',12); ylabel('r with reference','FontSize',12);
title('(A) Prior variance scaling','FontSize',12);
ylim([0.8 1.05]); grid on; box off; xline(1,'k--','LineWidth',1);

subplot(1,4,2);
plot(lambda_levels, sign_cov*100, 'bs-','MarkerFaceColor','b','LineWidth',2,'MarkerSize',8);
xlabel('\lambda','FontSize',12); ylabel('Sign consistency (%)','FontSize',12);
title('(B) Sign (\lambda)','FontSize',12);
ylim([0 105]); grid on; box off; xline(1,'k--','LineWidth',1);

subplot(1,4,3);
r_mean = nanmean(corr_mean_all,2); r_std = nanstd(corr_mean_all,0,2);
errorbar(alpha_levels_all, r_mean, r_std, 'ro-','MarkerFaceColor','r','LineWidth',2,'MarkerSize',8,'CapSize',6);
xlabel('\alpha','FontSize',12); ylabel('r with reference','FontSize',12);
title('(C) Prior mean perturbation','FontSize',12);
ylim([0 1.05]); grid on; box off;

subplot(1,4,4);
s_mean = nanmean(sign_mean_all,2)*100; s_std = nanstd(sign_mean_all,0,2)*100;
errorbar(alpha_levels_all, s_mean, s_std, 'ms-','MarkerFaceColor','m','LineWidth',2,'MarkerSize',8,'CapSize',6);
xlabel('\alpha','FontSize',12); ylabel('Sign consistency (%)','FontSize',12);
title('(D) Sign (\alpha)','FontSize',12);
ylim([0 105]); grid on; box off;

sgtitle(sprintf('S16 RegLin-Mean 12n prior sensitivity (N_{mc}=%d, T=672)', N_mc), ...
    'FontSize',13,'FontWeight','bold');
exportgraphics(fig, out_fig, 'Resolution',300);
fprintf('Figure saved: %s\nDone.\n', out_fig);
close(fig);
