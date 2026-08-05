%% run_splithalf_iter2_halfB_reglin_mean12n.m
%% Split-half iter2 — HALF B only: S13,S15,S17 (n=3)
%% Designed to run on server 86 in parallel with Half A on server 89.
%% 2026 - IMAG-26-0111

clear; clc;
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
out_mat  = fullfile(data_dir,'splithalf_iter2_halfB_reglin_mean12n.mat');
n = 12; Np = n*n;
mask_idx = find(~eye(n));
M_peb.Q  = 'single';

subjects = [13 15 17];
fprintf('\n=== Split-half iter2: Half B (S%s) on server 86 ===\n', num2str(subjects));
[BMA_B, PEB_B_i2, Ep_B, Pp_B] = run_iter2_half(subjects, data_dir, n, Np, M_peb);

Pp_B_vec = Pp_B(mask_idx);
Ep_B_vec = Ep_B(mask_idx);
sigB     = Pp_B_vec > 0.975;

fprintf('\n=== Half-B Results ===\n');
fprintf('Pp>0.975: %d | Pp>0.95: %d\n', sum(sigB), sum(Pp_B_vec>0.95));
fprintf('PEB_B_i2 F=%.2f\n', PEB_B_i2.F);

sh2B.BMA_B   = BMA_B;
sh2B.PEB_B   = PEB_B_i2;
sh2B.Ep_B    = Ep_B;
sh2B.Pp_B    = Pp_B;
sh2B.subjects = subjects;
save(out_mat, 'sh2B', '-v7.3');
fprintf('Saved: %s\n', out_mat);
fprintf('=== DONE ===\n');

function [BMA, PEB_i2, Ep_mat, Pp_mat] = run_iter2_half(subjects, data_dir, n, Np, M_peb)
    fprintf('  Loading %d original DCMs...\n', numel(subjects));
    GCM0 = cell(numel(subjects),1);
    for si = 1:numel(subjects)
        s = subjects(si);
        f = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', s));
        tmp = load(f,'DCM_est');
        GCM0{si} = tmp.DCM_est;
    end

    fprintf('  First-pass PEB (n=%d)...\n', numel(subjects));
    [PEB1, ~] = spm_dcm_peb(GCM0, M_peb, {'A'});
    fprintf('  PEB1 F=%.2f\n', PEB1.F);

    A_group = reshape(full(real(PEB1.Ep(1:Np,1))), n, n);
    if isfield(PEB1,'Ce') && ~isempty(PEB1.Ce)
        Ce_diag  = diag(full(real(PEB1.Ce)));
        pC_A_vec = max(Ce_diag(1:Np), 1e-4);
    else
        A_stack = zeros(n, n, numel(GCM0));
        for i = 1:numel(GCM0)
            A_stack(:,:,i) = GCM0{i}.Ep.A;
        end
        pC_A_vec = max(reshape(var(A_stack,0,3), Np, 1), 1e-4);
    end
    pC_A_mat = reshape(pC_A_vec, n, n);

    GCM_i2 = cell(numel(subjects),1);
    for si = 1:numel(subjects)
        s = subjects(si);
        fprintf('  Reinverting S%d...', s);
        sig_file = fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', s));
        d2 = load(sig_file,'Y_mean');
        Y  = detrend(d2.Y_mean);
        if size(Y,2) ~= n, Y = Y'; end

        DCM_orig = GCM0{si};
        pE = DCM_orig.M.pE;  pE.A = A_group;
        pC = DCM_orig.M.pC;  pC.A = pC_A_mat;

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
            DCM_est    = spm_dcm_calcium_csd(DCM2);
            GCM_i2{si} = DCM_est;
            fprintf(' F=%.2f\n', DCM_est.F);
        catch ME
            fprintf(' FAILED: %s (keeping first-pass)\n', ME.message);
            GCM_i2{si} = GCM0{si};
        end
    end

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
