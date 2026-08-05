% Diagnostic: why did S14's r_iter2 (0.0899) come out lower than r_flat (0.1238)?
% Read-only: loads already-cached h1/h2 DCMs, recomputes the mini-PEB (cheap),
% and inspects what changed between raw and reinverted A-matrices.
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
n = 12;
mask_idx = find(~eye(n));
s = 14;

d1 = load(fullfile(data_dir, sprintf('splithalf_within_reglin_mean_12n_s%d_h1.mat', s)), 'DCM_est');
d2 = load(fullfile(data_dir, sprintf('splithalf_within_reglin_mean_12n_s%d_h2.mat', s)), 'DCM_est');
DCM_half = {d1.DCM_est, d2.DCM_est};

Ep_h1 = full(DCM_half{1}.Ep.A);
Ep_h2 = full(DCM_half{2}.Ep.A);
r_flat = corr(Ep_h1(mask_idx), Ep_h2(mask_idx));
fprintf('r_flat (raw, recomputed) = %.4f\n', r_flat);

% Posterior uncertainty of each raw half (off-diag A block only)
Cp1 = full(diag(DCM_half{1}.Cp)); Cp2 = full(diag(DCM_half{2}.Cp));
sd1 = sqrt(Cp1(1:n*n)); sd1 = sd1(mask_idx);
sd2 = sqrt(Cp2(1:n*n)); sd2 = sd2(mask_idx);
fprintf('Raw posterior SD (off-diag A): half1 mean=%.4f median=%.4f | half2 mean=%.4f median=%.4f\n', ...
    mean(sd1), median(sd1), mean(sd2), median(sd2));
fprintf('Raw F: half1=%.2f  half2=%.2f  (delta=%.2f)\n', DCM_half{1}.F, DCM_half{2}.F, DCM_half{1}.F-DCM_half{2}.F);

rng(0);
M_peb.Q  = 'single';
[PEB_s, RCM_s] = spm_dcm_peb(DCM_half, M_peb, {'A'});
fprintf('class(PEB_s)=%s  numel(PEB_s)=%d\n', class(PEB_s), numel(PEB_s));
for k = 1:numel(PEB_s)
    fprintf('  PEB_s(%d).F=%.2f\n', k, PEB_s(k).F);
end

Ep_r1 = full(RCM_s{1}.Ep.A);
Ep_r2 = full(RCM_s{2}.Ep.A);
r_iter2 = corr(Ep_r1(mask_idx), Ep_r2(mask_idx));
fprintf('r_iter2 (reinverted, recomputed) = %.4f\n', r_iter2);

% How far did each half move from its raw estimate toward the shared prior?
shift1 = Ep_r1(mask_idx) - Ep_h1(mask_idx);
shift2 = Ep_r2(mask_idx) - Ep_h2(mask_idx);
fprintf('Shift magnitude (raw->reinverted): half1 mean|shift|=%.4f  half2 mean|shift|=%.4f\n', ...
    mean(abs(shift1)), mean(abs(shift2)));
fprintf('Correlation between the two shift vectors: %.4f (positive=pulled same direction, negative=pulled apart)\n', ...
    corr(shift1, shift2));

% Sign agreement before vs after, on the top-10 |raw mean(h1,h2)| edges
mean_raw = (Ep_h1(mask_idx) + Ep_h2(mask_idx)) / 2;
[~, ord] = sort(abs(mean_raw), 'descend');
top10 = ord(1:10);
fprintf('\nTop-10 largest-magnitude edges (raw mean): raw_h1 raw_h2 | reinv_h1 reinv_h2\n');
h1v = Ep_h1(mask_idx); h2v = Ep_h2(mask_idx);
r1v = Ep_r1(mask_idx); r2v = Ep_r2(mask_idx);
for k = 1:10
    idx = top10(k);
    fprintf('  edge %4d:  raw %+.3f %+.3f  ->  reinv %+.3f %+.3f  %s\n', ...
        idx, h1v(idx), h2v(idx), r1v(idx), r2v(idx), ...
        ternary(sign(h1v(idx))==sign(h2v(idx)) && sign(r1v(idx))~=sign(r2v(idx)), '<-- sign agreement LOST', ...
        ternary(sign(h1v(idx))~=sign(h2v(idx)) && sign(r1v(idx))==sign(r2v(idx)), '<-- sign agreement GAINED','')));
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
