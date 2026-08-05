% 4-way BMS: {RegHill, RegLin} x {PC1, Mean}, all M1-flat.
% Requires all 28 files (7 subjects x 4 models).
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
subjects  = 12:18;
n_subj    = numel(subjects);

model_tags   = {'sc_flat_reghill_12n',      'sc_flat_reglin_12n', ...
                'sc_flat_reghill_mean12n',   'sc_flat_reglin_mean12n'};
model_labels = {'RegHill-PC1','RegLin-PC1','RegHill-mean','RegLin-mean'};

fprintf('\n=== 4-way BMS: obs model x signal type (12n, M1-flat) ===\n');
F = zeros(n_subj, 4);
for mi = 1:4
    for si = 1:n_subj
        s  = subjects(si);
        fp = fullfile(data_dir, sprintf('subject_%d_DCM_%s.mat', s, model_tags{mi}));
        if ~exist(fp,'file')
            error('Missing: %s', fp);
        end
        tmp = load(fp,'DCM_est');
        F(si,mi) = tmp.DCM_est.F;
        fprintf('  %s S%d: F=%.2f\n', model_labels{mi}, s, F(si,mi));
    end
end

fprintf('\n--- Free energy table [subjects x models] ---\n');
fprintf('%-6s', 'Subj');
for mi = 1:4; fprintf('  %16s', model_labels{mi}); end; fprintf('\n');
for si = 1:n_subj
    fprintf('S%-5d', subjects(si));
    for mi = 1:4; fprintf('  %16.2f', F(si,mi)); end; fprintf('\n');
end

fprintf('\n--- 4-way random-effects BMS ---\n');
[alpha_r, exp_r, xp_r, pxp_r] = spm_BMS(F);
for mi = 1:4
    fprintf('  %s: alpha=%.2f  r=%.3f  pxp=%.3f\n', ...
        model_labels{mi}, alpha_r(mi), exp_r(mi), pxp_r(mi));
end
[~,best] = max(pxp_r);
fprintf('Winner: %s (PXP=%.3f)\n', model_labels{best}, pxp_r(best));

fprintf('\n--- Obs model BMS (pool signal types) ---\n');
F_obs = [max(F(:,[1,3]),[],2), max(F(:,[2,4]),[],2)];
[~,~,~,pxp_obs] = spm_BMS(F_obs);
fprintf('  RegHill (best signal): pxp=%.3f\n', pxp_obs(1));
fprintf('  RegLin  (best signal): pxp=%.3f\n', pxp_obs(2));

fprintf('\n--- Signal BMS (pool obs models) ---\n');
F_sig = [max(F(:,[1,2]),[],2), max(F(:,[3,4]),[],2)];
[~,~,~,pxp_sig] = spm_BMS(F_sig);
fprintf('  PC1  (best obs): pxp=%.3f\n', pxp_sig(1));
fprintf('  Mean (best obs): pxp=%.3f\n', pxp_sig(2));

fprintf('\n--- dF vs original global-linear M1 ---\n');
F_orig = nan(n_subj,1);
for si = 1:n_subj
    s = subjects(si);
    fp = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_12n.mat', s));
    if exist(fp,'file')
        tmp = load(fp,'DCM_est');
        F_orig(si) = tmp.DCM_est.F;
    end
end
if ~any(isnan(F_orig))
    fprintf('%-6s  %12s', 'Subj', 'GlobLin');
    for mi = 1:4; fprintf('  %16s', ['d(' model_labels{mi} ')']); end; fprintf('\n');
    for si = 1:n_subj
        fprintf('S%-5d  %12.1f', subjects(si), F_orig(si));
        for mi = 1:4; fprintf('  %+16.1f', F(si,mi)-F_orig(si)); end; fprintf('\n');
    end
    fprintf('Mean dF:  %12s', '');
    for mi = 1:4; fprintf('  %+16.1f', mean(F(:,mi)-F_orig)); end; fprintf('\n');
end

save(fullfile(data_dir,'bms_4way_12n.mat'), ...
    'F','alpha_r','exp_r','xp_r','pxp_r','model_labels','subjects', ...
    'F_obs','F_sig','pxp_obs','pxp_sig','-v7.3');
fprintf('\nSaved bms_4way_12n.mat\nDone.\n');
