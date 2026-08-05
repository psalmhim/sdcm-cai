%% plot_mc_12node_fig.m
% Regenerate FigS_mc_12node.png from saved results (no re-simulation).
% Run this after run_mc_simulation_12node.m has already been executed.

clear; clc;

home_dir = getenv('HOME');
spm_path = fullfile(home_dir, 'Dropbox/matlabwork/spm25');
dcm_path = fullfile(home_dir, 'Dropbox/matlabwork/mnet0.92/dcmcai');
addpath(spm_path); addpath(dcm_path);
cd(dcm_path);

in_mat  = './zebra/mc_12node_results.mat';
out_fig = './figures/FigS_mc_12node.png';

d = load(in_mat);
A_true     = d.A_true;
Aest_all   = d.Aest_all;
r_vals     = d.r_vals;
rmse_vals  = d.rmse_vals;
sign_acc   = d.sign_acc;
node_names = d.node_names;
N_mc       = d.N_mc;
n          = size(A_true, 1);

offdiag_mask = ~eye(n);
Amean        = nanmean(Aest_all, 3);

fprintf('Loaded: %s\n', in_mat);
fprintf('r = %.3f ± %.3f\n', nanmean(r_vals), nanstd(r_vals));
fprintf('RMSE = %.3f ± %.3f\n', nanmean(rmse_vals), nanstd(rmse_vals));
fprintf('Sign acc = %.1f ± %.1f%%\n', nanmean(sign_acc)*100, nanstd(sign_acc)*100);

%% Figure — doubled size; panel C gets extra width for 12×12 matrix + colorbar
fig = figure('Position', [100 100 2800 960], 'Color', 'w', 'Visible', 'off');

% Unequal column widths: A=0.25, B=0.25, C=0.40
left_margin = 0.07; right_margin = 0.04; gap = 0.06;
top = 0.13; bot = 0.13;
panel_h = 1 - top - bot;
w_A = 0.25; w_B = 0.25;
w_C = 1 - left_margin - right_margin - w_A - w_B - 2*gap;
x_A = left_margin;
x_B = x_A + w_A + gap;
x_C = x_B + w_B + gap;

% Panel A: scatter true vs estimated
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
title(sprintf('(A) 12-node recovery\nr = %.3f \\pm %.3f', ...
    nanmean(r_vals), nanstd(r_vals)), 'FontSize', 13);
axis square; grid on; box off; set(ax1, 'FontSize', 12);

% Panel B: per-run r bar chart
ax2 = axes('Position', [x_B bot w_B panel_h]);
bar(1:N_mc, r_vals, 0.7, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
hold on;
yline(nanmean(r_vals), 'r--', 'LineWidth', 2);
xlabel('Monte Carlo run', 'FontSize', 13);
ylabel('r (off-diagonal A)', 'FontSize', 13);
title(sprintf('(B) Per-run recovery\nmean r = %.3f, sign acc = %.1f%%', ...
    nanmean(r_vals), nanmean(sign_acc)*100), 'FontSize', 13);
ylim([-0.1 1]); grid on; box off; set(ax2, 'FontSize', 12);

% Panel C: ground-truth A matrix
ax3 = axes('Position', [x_C bot w_C panel_h]);
imagesc(A_true);
colormap(ax3, diverging_cmap(256));
cb = colorbar; cb.FontSize = 11; cb.Label.String = 'Connection strength';
clim([-max(abs(A_true(:))), max(abs(A_true(:)))]);
set(ax3, 'XTick', 1:n, 'XTickLabel', node_names, 'XTickLabelRotation', 45, ...
         'YTick', 1:n, 'YTickLabel', node_names, 'FontSize', 10);
title('(C) Ground truth A_{true} (12-node sparse)', 'FontSize', 13);
axis square; box off;

sgtitle(sprintf(['12-node sDCM-CaI: matched-model recovery  ' ...
    '(N_{mc}=%d, ~15%% density, sign acc=%.1f%%, RMSE=%.3f)'], ...
    N_mc, nanmean(sign_acc)*100, nanmean(rmse_vals)), ...
    'FontSize', 13, 'FontWeight', 'bold');

out_dir = fileparts(out_fig);
if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
exportgraphics(fig, out_fig, 'Resolution', 300);
fprintf('\nFigure saved: %s\n', out_fig);
close(fig);

function cmap = diverging_cmap(n)
half = floor(n/2);
r = [linspace(0.1,1,half), linspace(1,0.8,n-half)];
g = [linspace(0.3,1,half), linspace(1,0.1,n-half)];
b = [linspace(0.8,1,half), linspace(1,0.1,n-half)];
cmap = [r(:), g(:), b(:)];
end
