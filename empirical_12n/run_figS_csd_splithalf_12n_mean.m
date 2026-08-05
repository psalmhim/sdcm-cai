%% run_figS_csd_splithalf_12n_mean.m
% Regenerate FigS_csd_split_half.png using RegLin-Mean 12n data.
% CSD from meantrace_12n (mean signals); A-matrix from splithalf_within_reglin_mean_12n_results.mat.
% Replaces the old PC1-active 14n version.

maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
data_dir = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
out_fig  = fullfile(home_dir,'Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/FigS_csd_split_half.png');

subjects = 12:18;
n_subj   = numel(subjects);
n_nodes  = 12;
TR       = 0.5;
fs       = 1/TR;

%% Load A-matrix split-half results
sh = load(fullfile(data_dir,'splithalf_within_reglin_mean_12n_results.mat'));
r_A = sh.r_flat;   % 1x7, per-subject off-diagonal A-matrix split-half r

%% Compute CSD split-half from mean signals
r_auto = nan(1,n_subj);
r_cross = nan(1,n_subj);
T_half_all = nan(1,n_subj);

win_sec = 60;   % Welch window: 60 s
win_pts = round(win_sec * fs);
noverlap = round(win_pts * 0.5);
nfft = 2^nextpow2(win_pts * 2);

for si = 1:n_subj
    subj = subjects(si);
    fname = fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', subj));
    if ~exist(fname,'file')
        fprintf('Missing: %s\n', fname);
        continue;
    end
    d = load(fname,'Y_mean');
    sig = d.Y_mean;
    if size(sig,2) ~= n_nodes, sig = sig'; end
    T = size(sig,1);
    T_half = floor(T/2);
    T_half_all(si) = T_half;

    sig = detrend(sig);
    h1 = sig(1:T_half, :);
    h2 = sig(T_half+1:2*T_half, :);
    h1 = (h1 - mean(h1)) ./ (std(h1) + eps);
    h2 = (h2 - mean(h2)) ./ (std(h2) + eps);

    % Compute auto-spectra and cross-spectra for each half
    n_freq = nfft/2 + 1;
    auto1 = zeros(n_freq, n_nodes);
    auto2 = zeros(n_freq, n_nodes);
    cross1 = zeros(n_freq, n_nodes*(n_nodes-1)/2);
    cross2 = zeros(n_freq, n_nodes*(n_nodes-1)/2);

    for k = 1:n_nodes
        [p1, ~] = pwelch(h1(:,k), win_pts, noverlap, nfft, fs);
        [p2, ~] = pwelch(h2(:,k), win_pts, noverlap, nfft, fs);
        auto1(:,k) = p1;
        auto2(:,k) = p2;
    end

    idx = 0;
    for k = 1:n_nodes
        for l = k+1:n_nodes
            idx = idx+1;
            [c1, ~] = cpsd(h1(:,k), h1(:,l), win_pts, noverlap, nfft, fs);
            [c2, ~] = cpsd(h2(:,k), h2(:,l), win_pts, noverlap, nfft, fs);
            cross1(:,idx) = abs(c1);
            cross2(:,idx) = abs(c2);
        end
    end

    % Restrict to DCM band 0-0.15 Hz
    freqs = (0:nfft/2) * fs / nfft;
    band  = freqs <= 0.15;

    auto1_v  = auto1(band,:);
    auto2_v  = auto2(band,:);
    cross1_v = cross1(band,:);
    cross2_v = cross2(band,:);

    r_auto(si)  = corr(auto1_v(:),  auto2_v(:));
    r_cross(si) = corr(cross1_v(:), cross2_v(:));

    fprintf('S%d (T_half=%d): r_auto=%.3f  r_cross=%.3f  r_A=%.3f\n', ...
        subj, T_half, r_auto(si), r_cross(si), r_A(si));
end

fprintf('\nMean CSD auto=%.3f+/-%.3f  cross=%.3f+/-%.3f  A=%.3f+/-%.3f\n', ...
    nanmean(r_auto), nanstd(r_auto), nanmean(r_cross), nanstd(r_cross), ...
    nanmean(r_A), nanstd(r_A));

%% Figure
fig = figure('Position',[100 100 2400 800],'Color','w','Visible','off');

% Panel A: CSD cross-spectra vs A-matrix
subplot(1,3,1);
scatter(r_cross, r_A, 80, [0.2 0.4 0.8], 'filled'); hold on;
xl = [0 1]; plot(xl, xl, 'k--','LineWidth',1);
for si = 1:n_subj
    text(r_cross(si)+0.02, r_A(si), sprintf('S%d',subjects(si)), 'FontSize',9);
end
r_corr = corr(r_cross(:), r_A(:),'rows','complete');
text(0.05, 0.92, sprintf('r = %.3f', r_corr), 'FontSize',10, 'Units','data',...
    'FontSize',9);
xlabel('CSD cross-spectra split-half r','FontSize',11);
ylabel('A off-diagonal split-half r','FontSize',11);
title('(A) CSD vs A reliability','FontSize',11);
xlim([0 1]); ylim([-0.6 1]); grid on; box off;

% Panel B: per-subject bars
subplot(1,3,2);
x = 1:n_subj;
bw = 0.25;
bar(x - bw, r_auto,  bw*1.8, 'FaceColor',[0.6 0.6 0.6],'EdgeColor','none'); hold on;
bar(x,       r_cross, bw*1.8, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none');
bar(x + bw, r_A,     bw*1.8, 'FaceColor',[0.2 0.7 0.3],'EdgeColor','none');
yline(0,'k-','LineWidth',0.5);
set(gca,'XTick',1:n_subj,'XTickLabel',arrayfun(@(s)sprintf('S%d',s),subjects,'UniformOutput',false));
ylabel('Split-half r','FontSize',11);
title('(B) CSD and A reliability per subject','FontSize',11);
ylim([-0.6 1.05]); grid on; box off;
legend({'CSD auto-spectra','CSD cross-spectra','A off-diagonal'},'Location','southwest','FontSize',8);

% Panel C: A vs CSD colored by recording length
subplot(1,3,3);
scatter(r_cross, r_A, 80, T_half_all, 'filled'); hold on;
plot([0 1],[0 1],'k--','LineWidth',1);
for si = 1:n_subj
    text(r_cross(si)+0.02, r_A(si), sprintf('S%d',subjects(si)),'FontSize',9);
end
xlabel('CSD cross-spectra split-half r','FontSize',11);
ylabel('A off-diagonal split-half r','FontSize',11);
title('(C) Reliability vs recording length','FontSize',11);
cb = colorbar; cb.Label.String = 'Half recording (tp)'; colormap(gca,'cool');
xlim([0 1]); ylim([-0.6 1]); grid on; box off;

sgtitle('CSD-level vs A-matrix split-half reliability (RegLin-Mean 12n, mean signal)', ...
    'FontSize',12,'FontWeight','bold');

exportgraphics(fig, out_fig, 'Resolution',300);
fprintf('Saved: %s\nDone.\n', out_fig);
close(fig);
