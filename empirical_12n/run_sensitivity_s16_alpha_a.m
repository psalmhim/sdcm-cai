%% run_sensitivity_s16_alpha_a.m
% Prior mean perturbation, alpha=[0.05, 0.10], N_mc=10 each.
% Saves zebra/sensitivity_12n_s16_alpha_a.mat
% Run on server 89 in parallel with run_sensitivity_s16_alpha_b.m (alpha=[0.25,0.50] on 92).

maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir  = './zebra/';
out_mat   = fullfile(data_dir,'sensitivity_12n_s16_alpha_a.mat');

n           = 12;
N_mc        = 10;
alpha_levels = [0.05, 0.10];
TR = 0.5;

%% Clone setup from run_sc_flat_reglin_mean12n.m
d = load(fullfile(data_dir,'subject_16_meantrace_12n.mat'),'Y_mean');
signals = d.Y_mean;
if size(signals,2) ~= n; signals = signals'; end
Y = detrend(signals);

A = ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
options.nonlinear=0; options.two_state=0; options.stochastic=0;
options.induced=1; options.centre=1; options.modality='Ca';
options.Fdcm=[0.01 0.2]; options.Nmax=32; options.maxit=256;
options.precision=log(16); options.verbose=0; options.linear_obs=0;
options.custom_g = @spm_gx_calcium_linear;

[pE_ref,pC_ref,x0] = spm_dcm_calcium_priors(A,B,C,D,options);
pE_ref.alpha  = zeros(n,1); pC_ref.alpha  = ones(n,1) * 1/16;
pE_ref.beta_y = zeros(n,1); pC_ref.beta_y = ones(n,1) * 1/64;
pE_ref.Kd = 0; pC_ref.Kd = 1/128;
pE_ref.n  = 0; pC_ref.n  = 1/128;

pC_A = full(pC_ref.A);
offdiag = ~eye(n);

ref   = load(fullfile(data_dir,'subject_16_DCM_sc_flat_reglin_mean12n.mat'),'DCM_est');
Aref  = full(ref.DCM_est.Ep.A);
F_ref = ref.DCM_est.F;
fprintf('Reference F=%.4f  T=%d\n', F_ref, size(ref.DCM_est.Y.y,1));

DCM_tpl = struct();
DCM_tpl.a=A; DCM_tpl.b=B; DCM_tpl.c=C; DCM_tpl.d=D;
DCM_tpl.Y.y=Y; DCM_tpl.Y.dt=TR; DCM_tpl.Y.Q=[];
DCM_tpl.U.u=[]; DCM_tpl.U.dt=TR; DCM_tpl.M.nograph=1;
DCM_tpl.M.x=x0; DCM_tpl.options=options; DCM_tpl.options.maxnodes=max(16,n);

%% Alpha mean perturbation
fprintf('\n=== Alpha mean perturbation (alpha=[0.05,0.10], N_mc=%d) ===\n', N_mc);
corr_mean_all = nan(numel(alpha_levels), N_mc);
sign_mean_all = nan(numel(alpha_levels), N_mc);
F_mean_all    = nan(numel(alpha_levels), N_mc);

for ai = 1:numel(alpha_levels)
    alp = alpha_levels(ai);
    fprintf('  alpha=%.2f :', alp);
    for mc = 1:N_mc
        rng(ai * 1000 + mc);
        R = sign(randn(n,n)); R(~~eye(n)) = 0;
        pE_pert = pE_ref;
        pE_pert.A = pE_ref.A + alp * R .* sqrt(pC_A);
        DCM_run = DCM_tpl;
        DCM_run.options.pE = pE_pert;
        DCM_run.options.pC = pC_ref;
        t0 = tic;
        try
            est  = spm_dcm_calcium_csd(DCM_run);
            Aest = full(est.Ep.A);
            Fest = est.F;
        catch ME
            fprintf(' FAILED:%s', ME.message); Aest = nan(n,n); Fest = nan;
        end
        F_mean_all(ai,mc) = Fest;
        if ~any(isnan(Aest(:)))
            corr_mean_all(ai,mc) = corr(Aref(offdiag), Aest(offdiag));
            sign_mean_all(ai,mc) = mean(sign(Aref(offdiag)) == sign(Aest(offdiag)));
        end
        fprintf(' %.3f(%.0fm)', corr_mean_all(ai,mc), toc(t0)/60);
    end
    fprintf('\n  -> r=%.3f+/-%.3f  sign=%.1f+/-%.1f%%\n', ...
        nanmean(corr_mean_all(ai,:)), nanstd(corr_mean_all(ai,:)), ...
        nanmean(sign_mean_all(ai,:))*100, nanstd(sign_mean_all(ai,:))*100);
end

save(out_mat, 'alpha_levels','N_mc','n','corr_mean_all','sign_mean_all','F_mean_all','-v7.3');
fprintf('Saved: %s\nDone.\n', out_mat);
