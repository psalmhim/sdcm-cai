% Within-subject temporal split-half reliability for 12-node RegLin-Mean (flat SC, linear obs).
% Signals extracted from full-recording DCM (Y.y). No group data used.
% Design: for each subject, split recording in half -> fit DCM_h1 and DCM_h2
%         -> mini-PEB on the two halves -> reinverted h1 vs h2 (iter2).
% Results: r_flat and r_iter2 per subject, saved to splithalf_within_reglin_mean_12n_results.mat
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
subjects  = [12 13 14 15 16 17 18];
n = 12; TR = 0.5;

A = 0.001 * ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
options.nonlinear=0; options.two_state=0; options.stochastic=0;
options.induced=1; options.centre=1; options.modality='Ca';
options.Fdcm=[0.01 0.2]; options.Nmax=32; options.maxit=256;
options.precision=log(16); options.verbose=0; options.linear_obs=1;
options.custom_g = @spm_gx_calcium_linear;  % fitted alpha*Ca+beta (matches all other RegLin-Mean scripts)

[pE,pC,x0] = spm_dcm_calcium_priors(A,B,C,D,options);
pE.A = pE.A * 0;

% Per-region alpha/beta_y (matches run_sc_flat_reglin_mean12n.m exactly)
pE.alpha  = zeros(n,1);
pC.alpha  = ones(n,1) * 1/16;
pE.beta_y = zeros(n,1);
pC.beta_y = ones(n,1) * 1/64;
pE.Kd = 0; pC.Kd = 1/128;
pE.n  = 0; pC.n  = 1/128;

options.pE = pE; options.pC = pC;  % REQUIRED: spm_dcm_calcium_csd only honors overrides via options.pE/pC,
                                    % it recomputes its own pE/pC from scratch otherwise (ignores DCM.M.pE/pC)

M_peb.Q  = 'single';
mask_idx = find(~eye(n));

r_flat  = nan(1, numel(subjects));
r_iter2 = nan(1, numel(subjects));
failed  = false(1, numel(subjects));

for si = 1:numel(subjects)
    s = subjects(si);
    fprintf('\n=== Subject %d ===\n', s);
    try
        d = load(fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', s)), ...
                 'DCM_est');
        signals  = d.DCM_est.Y.y;  % T x 12
        T_full   = size(signals, 1);
        T_half   = floor(T_full / 2);
        fprintf('  T=%d -> T_half=%d\n', T_full, T_half);

        halves   = {signals(1:T_half, :), signals(T_half+1:2*T_half, :)};
        DCM_half = cell(2,1);

        for h = 1:2
            out = fullfile(data_dir, sprintf('splithalf_within_reglin_mean_12n_s%d_h%d.mat', s, h));
            Y_input = detrend(halves{h}) - mean(detrend(halves{h}));  % T_half x 12
            if exist(out, 'file')
                tmp = load(out, 'DCM_est');
                % Verify this is a RegLin-Mean inversion (not PC1 or Hill)
                if ~isfield(tmp.DCM_est,'options') || tmp.DCM_est.options.linear_obs ~= 1
                    error('OPTION MISMATCH: cached %s is not a linear_obs=1 (RegLin-Mean) inversion', out);
                end
                fprintf('  Half %d: cached [verified linear_obs=1]  ', h);
                DCM_half{h} = tmp.DCM_est;
                fprintf('F=%.2f\n', DCM_half{h}.F);
                continue;
            end
            fprintf('  Half %d: inverting...  ', h);
            Y = Y_input;

            DCM = struct();
            DCM.a=A; DCM.b=B; DCM.c=C; DCM.d=D;
            DCM.Y.y=Y; DCM.Y.dt=TR; DCM.Y.Q=[];
            DCM.U.u=[]; DCM.U.dt=TR; DCM.M.nograph=1;
            DCM.M.pE=pE; DCM.M.pC=pC; DCM.M.x=x0;
            DCM.M.n=size(x0,1); DCM.M.m=0; DCM.M.l=size(x0,2);
            DCM.options=options; DCM.options.maxnodes=max(16,n);
            DCM.name=sprintf('DCM_reglin_mean_12n_s%d_h%d', s, h);

            rng(0); t0 = tic;
            DCM_est = spm_dcm_calcium_csd(DCM);
            fprintf('F=%.2f  (%.1f min)\n', DCM_est.F, toc(t0)/60);
            DCM_half{h} = DCM_est;
            save(out, 'DCM_est', '-v7.3');
        end

        Ep_h1 = full(DCM_half{1}.Ep.A);
        Ep_h2 = full(DCM_half{2}.Ep.A);
        r_flat(si) = corr(Ep_h1(mask_idx), Ep_h2(mask_idx));
        fprintf('  r_flat: %.4f\n', r_flat(si));

        GCM_half = DCM_half;
        [PEB_s, RCM_s] = spm_dcm_peb(GCM_half, M_peb, {'A'});
        fprintf('  PEB_s F=%.2f\n', PEB_s.F);
        Ep_r1 = full(RCM_s{1}.Ep.A);
        Ep_r2 = full(RCM_s{2}.Ep.A);
        r_iter2(si) = corr(Ep_r1(mask_idx), Ep_r2(mask_idx));
        fprintf('  r_iter2: %.4f\n', r_iter2(si));
    catch ME
        failed(si) = true;
        fprintf('  *** SUBJECT %d FAILED: %s ***\n', s, ME.message);
        fprintf('  Skipping to next subject.\n');
        continue;
    end
end

fprintf('\n=== Split-half within-subject (12n RegLin-Mean) ===\n');
fprintf('Subjects: '); fprintf('S%d ', subjects); fprintf('\n');
fprintf('r_flat  : '); fprintf('%.3f ', r_flat);  fprintf('\n');
fprintf('r_iter2 : '); fprintf('%.3f ', r_iter2); fprintf('\n');
if any(failed)
    fprintf('FAILED subjects: '); fprintf('S%d ', subjects(failed)); fprintf('\n');
end
fprintf('r_flat  mean=%.3f  SD=%.3f  (n=%d of %d)\n', nanmean(r_flat),  nanstd(r_flat), sum(~isnan(r_flat)), numel(subjects));
fprintf('r_iter2 mean=%.3f  SD=%.3f  (n=%d of %d)\n', nanmean(r_iter2), nanstd(r_iter2), sum(~isnan(r_iter2)), numel(subjects));

save(fullfile(data_dir,'splithalf_within_reglin_mean_12n_results.mat'), ...
    'r_flat','r_iter2','subjects','failed','-v7.3');
fprintf('Saved splithalf_within_reglin_mean_12n_results.mat\nDone.\n');
