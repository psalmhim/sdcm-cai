% Within-subject temporal split-half iter2, PROPER full-reinversion version.
% Matches the established codebase pattern (run_peb_iter2_single_reglin_mean12n.m /
% run_splithalf_iter2_reglin_mean12n.m): first-pass PEB on the subject's own two
% halves -> extract empirical A_group/pC_A -> FULL nonlinear reinversion of each
% half via spm_dcm_calcium_csd with the empirical prior injected via BOTH
% DCM.M.pE/pC and DCM.options.pE/pC (the wrapper only honors options.pE/pC, but
% we set both defensively to match the established pattern exactly).
%
% This REPLACES the spm_dcm_peb closed-form RCM shortcut used in
% run_splithalf_within_reglin_mean_12n*.m for the iter2 step. The flat (raw)
% half-inversions are reused as-is from the already-cached .mat files -- only
% the iter2 reinversion is redone.
%
% Server-89 subset: subjects 12,14,16,18 (already has all flat halves cached).
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
subjects = [13 15 17];
n = 12; Np = n*n;
mask_idx = find(~eye(n));
M_peb.Q  = 'single';

r_flat  = nan(1, numel(subjects));
r_iter2 = nan(1, numel(subjects));
failed  = false(1, numel(subjects));

for si = 1:numel(subjects)
    s = subjects(si);
    fprintf('\n=== Subject %d (server 92, proper iter2) ===\n', s);
    try
        f1 = fullfile(data_dir, sprintf('splithalf_within_reglin_mean_12n_s%d_h1.mat', s));
        f2 = fullfile(data_dir, sprintf('splithalf_within_reglin_mean_12n_s%d_h2.mat', s));
        if ~exist(f1,'file') || ~exist(f2,'file')
            error('Missing cached flat half(s) for S%d', s);
        end
        d1 = load(f1, 'DCM_est'); d2 = load(f2, 'DCM_est');
        DCM_half = {d1.DCM_est, d2.DCM_est};

        % Reload raw signal and recompute the exact half-segments (matching how
        % the flat inversion built them), for a single clean detrend on reinversion.
        raw = load(fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', s)), 'Y_mean');
        signals = raw.Y_mean;
        if size(signals,2) ~= n, signals = signals'; end
        T_full = size(signals,1); T_half = floor(T_full/2);
        halves_raw = {signals(1:T_half,:), signals(T_half+1:2*T_half,:)};

        for h = 1:2
            if ~isfield(DCM_half{h},'options') || DCM_half{h}.options.linear_obs ~= 1
                error('OPTION MISMATCH: half %d is not a linear_obs=1 (RegLin-Mean) inversion', h);
            end
        end

        Ep_h1 = full(DCM_half{1}.Ep.A);
        Ep_h2 = full(DCM_half{2}.Ep.A);
        r_flat(si) = corr(Ep_h1(mask_idx), Ep_h2(mask_idx));
        fprintf('  r_flat: %.4f\n', r_flat(si));

        % Step 1-2: empirical prior directly from the two raw half posteriors.
        % NOTE: spm_dcm_peb(DCM_half, M_peb, {'A'}) with only 2 units returns
        % PEB1 as a 1x2 struct ARRAY (not a scalar PEB struct as it does for
        % n>=3), so PEB1.Ep(1:Np,1) throws a comma-separated-list indexing
        % error. Rather than depend on that n=2 edge-case behavior, compute
        % the group mean/variance directly -- this is exactly the same
        % fallback formula the established cross-subject scripts use when
        % PEB.Ce is unavailable, so it's already a validated approach.
        A_group  = (Ep_h1 + Ep_h2) / 2;
        A_stack  = zeros(n, n, 2);
        A_stack(:,:,1) = Ep_h1; A_stack(:,:,2) = Ep_h2;
        pC_A_vec = max(reshape(var(A_stack,0,3), Np, 1), 1e-4);
        pC_A_mat = reshape(pC_A_vec, n, n);

        % Step 3: full reinversion of each half under the empirical prior
        DCM_r = cell(2,1);
        for h = 1:2
            fprintf('  Reinverting half %d...', h);
            DCM_orig = DCM_half{h};
            pE = DCM_orig.M.pE;  pE.A = A_group;
            pC = DCM_orig.M.pC;  pC.A = pC_A_mat;

            Y_input = detrend(halves_raw{h}) - mean(detrend(halves_raw{h}));  % matches flat-inversion prep exactly

            DCM2 = struct();
            DCM2.a = DCM_orig.a; DCM2.b = DCM_orig.b; DCM2.c = DCM_orig.c; DCM2.d = DCM_orig.d;
            DCM2.Y.y = Y_input; DCM2.Y.dt = DCM_orig.Y.dt; DCM2.Y.Q = [];  % fresh single detrend, let wrapper rescale once
            DCM2.U.u = []; DCM2.U.dt = DCM_orig.U.dt;
            DCM2.M    = DCM_orig.M;  DCM2.M.pE = pE;  DCM2.M.pC = pC;
            DCM2.options = DCM_orig.options;
            DCM2.options.pE = pE;  DCM2.options.pC = pC;
            DCM2.options.maxnodes = max(16, n);
            DCM2.name = sprintf('splithalf_within_iter2proper_s%d_h%d', s, h);

            rng(h);
            DCM_est = spm_dcm_calcium_csd(DCM2);
            fprintf(' F=%.2f\n', DCM_est.F);
            DCM_r{h} = DCM_est;
            out = fullfile(data_dir, sprintf('splithalf_within_iter2proper_reglin_mean12n_s%d_h%d.mat', s, h));
            save(out, 'DCM_est', 'A_group', 'pC_A_vec', '-v7.3');
        end

        Ep_r1 = full(DCM_r{1}.Ep.A);
        Ep_r2 = full(DCM_r{2}.Ep.A);
        r_iter2(si) = corr(Ep_r1(mask_idx), Ep_r2(mask_idx));
        fprintf('  r_iter2 (proper reinversion): %.4f\n', r_iter2(si));
    catch ME
        failed(si) = true;
        fprintf('  *** SUBJECT %d FAILED: %s ***\n', s, ME.message);
        continue;
    end
end

fprintf('\n=== Split-half iter2 (PROPER reinversion), server-92 subset ===\n');
fprintf('Subjects: '); fprintf('S%d ', subjects); fprintf('\n');
fprintf('r_flat  : '); fprintf('%.3f ', r_flat);  fprintf('\n');
fprintf('r_iter2 : '); fprintf('%.3f ', r_iter2); fprintf('\n');
if any(failed)
    fprintf('FAILED subjects: '); fprintf('S%d ', subjects(failed)); fprintf('\n');
end
save(fullfile(data_dir,'splithalf_within_iter2proper_s92_results.mat'), ...
    'r_flat','r_iter2','subjects','failed','-v7.3');
fprintf('Saved splithalf_within_iter2proper_s92_results.mat\nDone.\n');
