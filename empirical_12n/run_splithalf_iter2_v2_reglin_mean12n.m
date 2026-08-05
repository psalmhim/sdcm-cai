
%% run_splithalf_iter2_v2_reglin_mean12n.m
%% Split-half iter2 — correct version: no final group PEB.
%% Each half: first-pass PEB -> empirical prior -> reinvert subjects -> mean Ep
%% Compare mean_A vs mean_B directly.
%% 2026 - IMAG-26-0111

clear; clc;
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
out_mat  = fullfile(data_dir,'splithalf_iter2_v2_reglin_mean12n.mat');
n = 12; Np = n*n;
mask_idx = find(~eye(n));
M_peb.Q  = 'single';

subj_A = [12 14 16 18];
subj_B = [13 15 17];

fprintf('\n=== Split-half iter2 v2 (no final group PEB) ===\n');

[mean_A, Ep_subj_A] = run_half(subj_A, data_dir, n, Np, M_peb, 'A');
[mean_B, Ep_subj_B] = run_half(subj_B, data_dir, n, Np, M_peb, 'B');

%% Metrics
a_vec = mean_A(mask_idx);
b_vec = mean_B(mask_idx);
r_all = corr(a_vec, b_vec);

sigA = sign(a_vec);
sigB = sign(b_vec);
sign_cons = mean(sigA == sigB);

fprintf('\n=== RESULTS ===\n');
fprintf('r_all (all 132 off-diag):  %.3f\n', r_all);
fprintf('sign_cons:                  %.3f\n', sign_cons);

%% Save
sh2v2.mean_A    = mean_A;
sh2v2.mean_B    = mean_B;
sh2v2.Ep_subj_A = Ep_subj_A;
sh2v2.Ep_subj_B = Ep_subj_B;
sh2v2.r_all     = r_all;
sh2v2.sign_cons = sign_cons;
sh2v2.subj_A    = subj_A;
sh2v2.subj_B    = subj_B;
save(out_mat, 'sh2v2', '-v7.3');
fprintf('Saved: %s\n', out_mat);
fprintf('=== DONE ===\n');

function [mean_Ep, Ep_all] = run_half(subjects, data_dir, n, Np, M_peb, label)
    fprintf('\n--- Half %s (S%s) ---\n', label, num2str(subjects));
    GCM0 = cell(numel(subjects),1);
    for si = 1:numel(subjects)
        s = subjects(si);
        tmp = load(fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat',s)),'DCM_est');
        GCM0{si} = tmp.DCM_est;
    end

    fprintf('  First-pass PEB...\n');
    [PEB1, ~] = spm_dcm_peb(GCM0, M_peb, {'A'});
    fprintf('  PEB F=%.2f\n', PEB1.F);

    A_group  = reshape(full(real(PEB1.Ep(1:Np,1))), n, n);
    Ce_diag  = diag(full(real(PEB1.Ce)));
    pC_A_vec = max(Ce_diag(1:Np), 1e-4);
    pC_A_mat = reshape(pC_A_vec, n, n);

    Ep_all = zeros(n, n, numel(subjects));
    for si = 1:numel(subjects)
        s = subjects(si);
        fprintf('  Reinverting S%d...', s);
        d2  = load(fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat',s)),'Y_mean');
        Y   = detrend(d2.Y_mean);
        if size(Y,2) ~= n, Y = Y'; end

        DCM_orig = GCM0{si};
        pE = DCM_orig.M.pE;  pE.A = A_group;
        pC = DCM_orig.M.pC;  pC.A = pC_A_mat;

        DCM2 = struct();
        DCM2.a=DCM_orig.a; DCM2.b=DCM_orig.b; DCM2.c=DCM_orig.c; DCM2.d=DCM_orig.d;
        DCM2.Y.y=Y; DCM2.Y.dt=DCM_orig.Y.dt; DCM2.Y.Q=[];
        DCM2.U.u=[]; DCM2.U.dt=DCM_orig.U.dt;
        DCM2.M=DCM_orig.M; DCM2.M.pE=pE; DCM2.M.pC=pC;
        DCM2.options=DCM_orig.options; DCM2.options.pE=pE; DCM2.options.pC=pC;
        DCM2.options.maxnodes=max(16,n);
        rng(si);
        try
            DCM_est = spm_dcm_calcium_csd(DCM2);
            Ep_all(:,:,si) = DCM_est.Ep.A;
            fprintf(' F=%.2f\n', DCM_est.F);
        catch ME
            fprintf(' FAILED: %s\n', ME.message);
            Ep_all(:,:,si) = GCM0{si}.Ep.A;
        end
    end
    mean_Ep = mean(Ep_all, 3);
    fprintf('  mean_Ep computed (no group PEB)\n');
end
