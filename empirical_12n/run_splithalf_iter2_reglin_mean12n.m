%% run_splithalf_iter2_reglin_mean12n.m
%% Fully independent split-half iter2 for RegLin-Mean 12n.
%% Each half: first-pass PEB -> subject reinversion -> iter2 PEB -> BMA
%%
%% Half A: S12,S14,S16,S18 (n=4)   Half B: S13,S15,S17 (n=3)

clear; clc;
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
out_mat  = fullfile(data_dir,'splithalf_iter2_reglin_mean12n.mat');
n = 12; Np = n*n;
mask_idx = find(~eye(n));
M_peb.Q  = 'single';

subj_A = [12 14 16 18];
subj_B = [13 15 17];

%% ── Run halves ────────────────────────────────────────────────────────────
fprintf('\n=== Split-half iter2: Half A (S%s) ===\n', num2str(subj_A));
[BMA_A, PEB_A_i2, Ep_A, Pp_A] = run_iter2_half(subj_A, data_dir, n, Np, M_peb, home_dir);

fprintf('\n=== Split-half iter2: Half B (S%s) ===\n', num2str(subj_B));
[BMA_B, PEB_B_i2, Ep_B, Pp_B] = run_iter2_half(subj_B, data_dir, n, Np, M_peb, home_dir);

%% ── Metrics ───────────────────────────────────────────────────────────────
Ep_A_vec = Ep_A(mask_idx);
Ep_B_vec = Ep_B(mask_idx);
Pp_A_vec = Pp_A(mask_idx);
Pp_B_vec = Pp_B(mask_idx);

r_all = corr(Ep_A_vec, Ep_B_vec);

% Full-sample reference Pp (iter2 PEB on all 7 subjects)
ref_file = fullfile(data_dir,'PEB_iter2_sc_flat_reglin_mean12n_results.mat');
if exist(ref_file,'file')
    ref = load(ref_file,'Pp_iter2');
    if isfield(ref,'Pp_iter2')
        Pp_full_off = ref.Pp_iter2(mask_idx);
        sig_idx_full = mask_idx(Pp_full_off > 0.975);
        if numel(sig_idx_full) > 1
            r_sig      = corr(Ep_A(sig_idx_full), Ep_B(sig_idx_full));
            n_sig_full = numel(sig_idx_full);
        else
            r_sig = NaN; n_sig_full = 0;
        end
    else
        r_sig = NaN; n_sig_full = NaN;
    end
else
    r_sig = NaN; n_sig_full = NaN;
    fprintf('Warning: full-sample iter2 PEB not found for r_sig\n');
end

sigA = Pp_A_vec > 0.975;
sigB = Pp_B_vec > 0.975;
inter   = double(sum(sigA & sigB));
union_n = double(sum(sigA | sigB));
jaccard = inter / max(union_n, 1);

shared = find(sigA & sigB);
if ~isempty(shared)
    same_sign = sum(sign(Ep_A_vec(shared)) == sign(Ep_B_vec(shared)));
    sign_cons = same_sign / numel(shared);
else
    same_sign = 0; sign_cons = NaN;
end

fprintf('\n=== RESULTS (iter2 split-half) ===\n');
fprintf('r (all 132 off-diag): %.3f\n',      r_all);
fprintf('r (full-sample sig):  %.3f (n=%d)\n', r_sig, n_sig_full);
fprintf('Jaccard:               %.3f (%d/%d)\n', jaccard, inter, union_n);
fprintf('Sign consistency:      %.3f (%d/%d)\n', sign_cons, same_sign, numel(shared));
fprintf('Half-A Pp>0.975: %d | Half-B Pp>0.975: %d\n', sum(sigA), sum(sigB));

%% ── Save ──────────────────────────────────────────────────────────────────
sh2.r_all      = r_all;
sh2.r_sig      = r_sig;
sh2.n_sig_full = n_sig_full;
sh2.jaccard    = jaccard;
sh2.sign_cons  = sign_cons;
sh2.inter      = inter;
sh2.union_n    = union_n;
sh2.Ep_A       = Ep_A;  sh2.Pp_A = Pp_A;
sh2.Ep_B       = Ep_B;  sh2.Pp_B = Pp_B;
sh2.PEB_A_i2   = PEB_A_i2;  sh2.BMA_A = BMA_A;
sh2.PEB_B_i2   = PEB_B_i2;  sh2.BMA_B = BMA_B;
sh2.subj_A     = subj_A;
sh2.subj_B     = subj_B;
save(out_mat, 'sh2', '-v7.3');
fprintf('\nSaved: %s\n', out_mat);
fprintf('=== DONE ===\n');

%% ════════════════════════════════════════════════════════════════════════════
%% LOCAL FUNCTION
%% ════════════════════════════════════════════════════════════════════════════
function [BMA, PEB_i2, Ep_mat, Pp_mat] = run_iter2_half(subjects, data_dir, n, Np, M_peb, home_dir)
    fprintf('  Loading %d original DCMs...\n', numel(subjects));
    GCM0 = cell(numel(subjects),1);
    for si = 1:numel(subjects)
        s = subjects(si);
        f = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', s));
        tmp = load(f,'DCM_est');
        GCM0{si} = tmp.DCM_est;
    end

    % Step 1: first-pass PEB for this half only
    fprintf('  First-pass PEB (n=%d)...\n', numel(subjects));
    [PEB1, ~] = spm_dcm_peb(GCM0, M_peb, {'A'});
    fprintf('  PEB1 F=%.2f\n', PEB1.F);

    % Step 2: empirical prior from half-group PEB
    A_group = reshape(full(real(PEB1.Ep(1:Np,1))), n, n);
    if isfield(PEB1,'Ce') && ~isempty(PEB1.Ce)
        Ce_diag    = diag(full(real(PEB1.Ce)));
        pC_A_vec   = max(Ce_diag(1:Np), 1e-4);
    else
        A_stack = zeros(n, n, numel(GCM0));
        for i = 1:numel(GCM0)
            A_stack(:,:,i) = GCM0{i}.Ep.A;
        end
        pC_A_vec = max(reshape(var(A_stack,0,3), Np, 1), 1e-4);
    end
    pC_A_mat = reshape(pC_A_vec, n, n);

    % Step 3: reinvert each subject with empirical prior (iter2)
    GCM_i2 = cell(numel(subjects),1);
    for si = 1:numel(subjects)
        s = subjects(si);
        fprintf('  Reinverting S%d...', s);

        sig_file = fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', s));
        d2 = load(sig_file,'Y_mean');
        Y  = detrend(d2.Y_mean);
        if size(Y,2) ~= n, Y = Y'; end

        DCM_orig = GCM0{si};
        pE       = DCM_orig.M.pE;  pE.A = A_group;
        pC       = DCM_orig.M.pC;  pC.A = pC_A_mat;

        DCM2 = struct();
        DCM2.a = DCM_orig.a;  DCM2.b = DCM_orig.b;
        DCM2.c = DCM_orig.c;  DCM2.d = DCM_orig.d;
        DCM2.Y.y  = Y;  DCM2.Y.dt = DCM_orig.Y.dt;  DCM2.Y.Q = [];
        DCM2.U.u  = [];  DCM2.U.dt = DCM_orig.U.dt;
        DCM2.M    = DCM_orig.M;  DCM2.M.pE = pE;  DCM2.M.pC = pC;
        DCM2.options = DCM_orig.options;
        DCM2.options.pE = pE;  DCM2.options.pC = pC;
        DCM2.options.maxnodes = max(16, n);
        rng(si);
        try
            DCM_est   = spm_dcm_calcium_csd(DCM2);
            GCM_i2{si} = DCM_est;
            fprintf(' F=%.2f\n', DCM_est.F);
        catch ME
            fprintf(' FAILED: %s  (keeping first-pass)\n', ME.message);
            GCM_i2{si} = GCM0{si};
        end
    end

    % Step 4: iter2 group PEB + BMA for this half
    fprintf('  Iter2 group PEB (n=%d)...\n', numel(subjects));
    [PEB_i2, ~] = spm_dcm_peb(GCM_i2, M_peb, {'A'});
    fprintf('  PEB_i2 F=%.2f\n', PEB_i2.F);
    BMA = spm_dcm_peb_bmc(PEB_i2);

    Ep_mat = reshape(full(real(BMA.Ep(1:Np))), n, n);
    Pp_mat = reshape(full(real(BMA.Pp(1:Np))), n, n);
    mask   = ~eye(n,'logical');
    fprintf('  Pp>0.975: %d | Pp>0.95: %d\n', ...
        sum(Pp_mat(mask) > 0.975), sum(Pp_mat(mask) > 0.95));
end
