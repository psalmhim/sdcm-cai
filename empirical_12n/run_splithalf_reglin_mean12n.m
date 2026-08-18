%% run_splithalf_reglin_mean12n.m
%%
%% Split-half cross-subject reliability for RegLin-Mean 12n (primary model).
%% Split: A = S12,S14,S16,S18 (n=4)  |  B = S13,S15,S17 (n=3)
%%
%% Metrics:
%%   r_all     - Pearson r on all 132 off-diagonal iter2 Ep.A entries
%%   r_sig     - Pearson r restricted to full-group BMA Pp>0.975 edges
%%   Jaccard   - overlap of Pp>0.975 edges from independent PEB halves
%%   sign_cons - sign consistency of shared significant edges
%%
%% OUTPUT: zebra/splithalf_reglin_mean12n.mat
%%

clear; clc;
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir  = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
out_mat   = fullfile(data_dir,'splithalf_reglin_mean12n.mat');
subjects  = [12 13 14 15 16 17 18];
n         = 12;
idx_A     = [1 3 5 7];
idx_B     = [2 4 6];
mask_idx  = find(~eye(n));
M         = struct(); M.Q = 'single';

fprintf('=== Split-half reliability: RegLin-Mean 12n ===\n');

%% Load all 7 DCMs
GCM = cell(numel(subjects), 1);
for si = 1:numel(subjects)
    s  = subjects(si);
    fp = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', s));
    fprintf('Loading %s...\n', fp);
    tmp    = load(fp,'DCM_est');
    GCM{si} = tmp.DCM_est;
end
fprintf('\nAll 7 DCMs loaded.\n');

%% Full PEB on all 7 -> group-regularised iter2 posteriors (RCM)
fprintf('\n--- Full PEB (all 7 subjects) ---\n');
[PEB_full, RCM] = spm_dcm_peb(GCM, M, {'A'});
BMA_full  = spm_dcm_peb_bmc(PEB_full);
Pp_full   = full(BMA_full.Pp(1:n*n));
n_sig_full = sum(Pp_full(mask_idx) > 0.975);
fprintf('F_PEB = %.2f  |  Pp>0.975 edges: %d\n', PEB_full.F, n_sig_full);

%% Mean iter2 Ep.A per split
sum_A = zeros(n,n);
for i = idx_A; sum_A = sum_A + full(RCM{i}.Ep.A); end
mean_A_mat = sum_A / numel(idx_A);

sum_B = zeros(n,n);
for i = idx_B; sum_B = sum_B + full(RCM{i}.Ep.A); end
mean_B_mat = sum_B / numel(idx_B);

%% Pearson r - all 132 off-diagonal
r_all = corr(mean_A_mat(mask_idx), mean_B_mat(mask_idx));

%% Pearson r - restricted to full-group significant edges
sig_idx = intersect(mask_idx, find(Pp_full > 0.975));
if numel(sig_idx) > 1
    r_sig = corr(mean_A_mat(sig_idx), mean_B_mat(sig_idx));
else
    r_sig = NaN;
end
fprintf('r (all 132):    %.3f\n', r_all);
fprintf('r (sig edges):  %.3f  (n=%d)\n', r_sig, numel(sig_idx));

%% Independent PEB on each half -> Jaccard + sign consistency
fprintf('\n--- PEB split A (n=%d) ---\n', numel(idx_A));
[PEB_A, ~] = spm_dcm_peb(GCM(idx_A), M, {'A'});
BMA_A = spm_dcm_peb_bmc(PEB_A);
Pp_A  = full(BMA_A.Pp(1:n*n));
Ep_A  = full(BMA_A.Ep(1:n*n));
fprintf('F_PEB_A = %.2f  |  Pp>0.975: %d\n', PEB_A.F, sum(Pp_A(mask_idx)>0.975));

fprintf('\n--- PEB split B (n=%d) ---\n', numel(idx_B));
[PEB_B, ~] = spm_dcm_peb(GCM(idx_B), M, {'A'});
BMA_B = spm_dcm_peb_bmc(PEB_B);
Pp_B  = full(BMA_B.Pp(1:n*n));
Ep_B  = full(BMA_B.Ep(1:n*n));
fprintf('F_PEB_B = %.2f  |  Pp>0.975: %d\n', PEB_B.F, sum(Pp_B(mask_idx)>0.975));

sigA = Pp_A(mask_idx) > 0.975;
sigB = Pp_B(mask_idx) > 0.975;
inter   = double(sum(sigA & sigB));
union_n = double(sum(sigA | sigB));
jaccard = inter / max(union_n, 1);

shared_idx = mask_idx(sigA & sigB);
if numel(shared_idx) > 0
    same_sign = sum(sign(Ep_A(shared_idx)) == sign(Ep_B(shared_idx)));
    sign_cons = same_sign / numel(shared_idx);
else
    same_sign = 0;
    sign_cons = NaN;
end

fprintf('\n=== RESULTS ===\n');
fprintf('r (all 132 off-diag): %.3f\n',   r_all);
fprintf('r (Pp>0.975 edges):   %.3f  (n=%d)\n', r_sig, numel(sig_idx));
fprintf('Jaccard:               %.3f  (%d inter / %d union)\n', jaccard, inter, union_n);
fprintf('Sign consistency:      %.3f  (%d/%d)\n', sign_cons, same_sign, numel(shared_idx));
fprintf('Full-group N_sig:      %d\n', n_sig_full);

%% Save
sh.r_all      = r_all;
sh.r_sig      = r_sig;
sh.jaccard    = jaccard;
sh.sign_cons  = sign_cons;
sh.n_sig_full = n_sig_full;
sh.inter      = inter;
sh.union_n    = union_n;
sh.PEB_full   = PEB_full;
sh.BMA_full   = BMA_full;
sh.PEB_A      = PEB_A;  sh.BMA_A = BMA_A;
sh.PEB_B      = PEB_B;  sh.BMA_B = BMA_B;
sh.mean_A     = mean_A_mat;
sh.mean_B     = mean_B_mat;
sh.subjects   = subjects;
sh.idx_A      = idx_A;
sh.idx_B      = idx_B;

save(out_mat, 'sh', '-v7.3');
fprintf('\nSaved: %s\n', out_mat);
fprintf('=== DONE ===\n');
