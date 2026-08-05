% BMS across SC architectures for 12-node 12n PC1 signals.
% M1: flat-p16   M2: binary-SC   M3: graded-SC
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
subjects = 12:18;
n_subj   = numel(subjects);

model_tags   = {'sc_flat', 'sc_binary', 'sc_graded'};
model_labels = {'M1-flat','M2-binary-SC','M3-graded-SC'};
n_models = numel(model_tags);

fprintf('\n=== BMS: SC architecture comparison (12-node 12n PC1) ===\n');

F = zeros(n_subj, n_models);
for mi = 1:n_models
    for si = 1:n_subj
        s = subjects(si);
        fname = sprintf('subject_%d_DCM_%s_12n.mat', s, model_tags{mi});
        fpath = fullfile(data_dir, fname);
        if ~exist(fpath,'file')
            error('Missing: %s', fpath);
        end
        tmp = load(fpath, 'DCM_est');
        F(si,mi) = tmp.DCM_est.F;
        fprintf('  %s S%d: F=%.2f\n', model_labels{mi}, s, F(si,mi));
    end
end

fprintf('\nFree energy matrix [subjects x models]:\n');
for mi=1:n_models; fprintf('  %-14s', model_labels{mi}); end; fprintf('\n');
for si=1:n_subj
    fprintf('  S%d:', subjects(si));
    for mi=1:n_models; fprintf('  %10.2f', F(si,mi)); end
    fprintf('\n');
end

fprintf('\n--- Random-effects BMS ---\n');
[alpha, exp_r, xp, pxp, bor] = spm_BMS(F);
for mi=1:n_models
    fprintf('  %s: alpha=%.2f  r=%.3f  xp=%.3f  pxp=%.3f\n', ...
        model_labels{mi}, alpha(mi), exp_r(mi), xp(mi), pxp(mi));
end
[~, best] = max(pxp);
fprintf('Winner: %s (PXP=%.3f)\n', model_labels{best}, pxp(best));

dF = F - F(:,1);
fprintf('\n--- Delta-F relative to M1-flat ---\n');
for mi=2:n_models
    fprintf('  %s: mean dF=%.2f  [%.2f to %.2f]\n', ...
        model_labels{mi}, mean(dF(:,mi)), min(dF(:,mi)), max(dF(:,mi)));
end

save(fullfile(data_dir,'bms_sc_12n.mat'), ...
    'F','alpha','exp_r','xp','pxp','bor','model_labels','subjects','-v7.3');
fprintf('\nSaved bms_sc_12n.mat\nDone.\n');
