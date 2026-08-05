% Split-half reliability for all 4 obs models (M1-flat, 12n).
% Split: A = S12,14,16,18 (pos 1,3,5,7)  |  B = S13,15,17 (pos 2,4,6)
% Metrics: Pearson r on mean iter2 Ep.A (132 off-diag) + Jaccard (Pp>0.975)
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
subjects  = [12 13 14 15 16 17 18];
n         = 12;
M         = struct(); M.Q = 'single';
idx_A     = [1 3 5 7];   % S12,14,16,18
idx_B     = [2 4 6];     % S13,15,17
mask_idx  = find(~eye(n));   % 132 off-diagonal indices

model_tags   = {'sc_flat_reghill_12n',    'sc_flat_reglin_12n', ...
                'sc_flat_reghill_mean12n', 'sc_flat_reglin_mean12n'};
model_labels = {'RegHill-PC1','RegLin-PC1','RegHill-mean','RegLin-mean'};

results = struct();

for mi = 1:4
    tag   = model_tags{mi};
    label = model_labels{mi};
    fprintf('\n=== Split-half: %s ===\n', label);

    GCM = cell(numel(subjects), 1);
    for si = 1:numel(subjects)
        s  = subjects(si);
        fp = fullfile(data_dir, sprintf('subject_%d_DCM_%s.mat', s, tag));
        tmp = load(fp,'DCM_est');
        GCM{si} = tmp.DCM_est;
    end

    % Full PEB -> RCM (group-regularised iter2 DCMs)
    fprintf('  Full PEB (all 7)...\n');
    [PEB_full, RCM] = spm_dcm_peb(GCM, M, {'A'});
    BMA_full = spm_dcm_peb_bmc(PEB_full);
    Pp_full  = full(BMA_full.Pp(1:n*n));
    n_sig_full = sum(Pp_full(mask_idx) > 0.975);
    fprintf('  F_PEB=%.2f  Pp>0.975: %d\n', PEB_full.F, n_sig_full);

    % Mean iter2 Ep.A per split
    sum_A = zeros(n,n);
    for i = idx_A; sum_A = sum_A + full(RCM{i}.Ep.A); end
    mean_A = sum_A / numel(idx_A);

    sum_B = zeros(n,n);
    for i = idx_B; sum_B = sum_B + full(RCM{i}.Ep.A); end
    mean_B = sum_B / numel(idx_B);

    % Pearson r (all 132 off-diagonal)
    r_all = corr(mean_A(mask_idx), mean_B(mask_idx));

    % Pearson r (full-BMA significant edges only)
    sig_idx = intersect(mask_idx, find(Pp_full > 0.975));
    if numel(sig_idx) > 1
        r_sig = corr(mean_A(sig_idx), mean_B(sig_idx));
    else
        r_sig = NaN;
    end

    % Jaccard: independent PEB on each half
    fprintf('  PEB half A (%d subs)...\n', numel(idx_A));
    BMA_A = spm_dcm_peb_bmc(spm_dcm_peb(GCM(idx_A), M, {'A'}));
    Pp_A  = full(BMA_A.Pp(1:n*n));

    fprintf('  PEB half B (%d subs)...\n', numel(idx_B));
    BMA_B = spm_dcm_peb_bmc(spm_dcm_peb(GCM(idx_B), M, {'A'}));
    Pp_B  = full(BMA_B.Pp(1:n*n));

    sigA = Pp_A(mask_idx) > 0.975;
    sigB = Pp_B(mask_idx) > 0.975;
    inter   = double(sum(sigA & sigB));
    union_n = double(sum(sigA | sigB));
    jaccard = inter / max(union_n, 1);

    % Sign consistency among shared edges
    shared_idx = mask_idx(sigA & sigB);
    Ep_A_full  = full(BMA_A.Ep(1:n*n));
    Ep_B_full  = full(BMA_B.Ep(1:n*n));
    if numel(shared_idx) > 0
        same_sign = sum(sign(Ep_A_full(shared_idx)) == sign(Ep_B_full(shared_idx)));
        sign_cons = same_sign / numel(shared_idx);
    else
        sign_cons = NaN;
    end

    fprintf('  r (all 132):     %.3f\n', r_all);
    fprintf('  r (sig edges):   %.3f  (n=%d)\n', r_sig, numel(sig_idx));
    fprintf('  Jaccard:         %.3f  (%d shared / %d union)\n', jaccard, inter, union_n);
    fprintf('  Sign cons:       %.3f  (%d/%d same sign)\n', sign_cons, ...
            double(sum(sign(Ep_A_full(shared_idx)) == sign(Ep_B_full(shared_idx)))), numel(shared_idx));

    results(mi).label     = label;
    results(mi).r_all     = r_all;
    results(mi).r_sig     = r_sig;
    results(mi).jaccard   = jaccard;
    results(mi).sign_cons = sign_cons;
    results(mi).n_sig     = n_sig_full;
    results(mi).inter     = inter;
    results(mi).union_n   = union_n;
end

fprintf('\n=== Summary ===\n');
fprintf('%-16s  %7s  %7s  %7s  %8s  %5s\n', 'Model','r_all','r_sig','Jaccard','SignCons','N_sig');
for mi = 1:4
    fprintf('%-16s  %7.3f  %7.3f  %7.3f  %8.3f  %5d\n', ...
        results(mi).label, results(mi).r_all, results(mi).r_sig, ...
        results(mi).jaccard, results(mi).sign_cons, results(mi).n_sig);
end

save(fullfile(data_dir,'splithalf_4models_12n.mat'), 'results', '-v7.3');
fprintf('\nSaved splithalf_4models_12n.mat\nDone.\n');
