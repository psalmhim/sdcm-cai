% Batch extraction: Mean signals (12-node, k-means TeO split) for all subjects.
% Run ONCE before run_mean12n.m. Sequential (no parfor).
%
% Overwrites any existing subject_N_meantrace_12n.mat files.
% After this completes, run run_mean12n.m to re-run DCMs.
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));

subjects = [12 13 14 15 16 17 18];
w_sp     = 0.5;   % equal spatial / temporal weight (same as PC1 extraction)

fprintf('=== Mean signal extraction (k-means TeO, no z-score) ===\n');
fprintf('w_sp=%.2f, subjects: %s\n', w_sp, mat2str(subjects));

for si = 1:numel(subjects)
    s = subjects(si);
    fprintf('\n--- Subject %d (%d/%d) ---\n', s, si, numel(subjects));
    extract_mean_signals_12n_from_raw(s, w_sp);
end

fprintf('\n=== All extractions complete. Now run run_mean12n.m ===\n');
