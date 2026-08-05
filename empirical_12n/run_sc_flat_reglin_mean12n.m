function run_sc_flat_reglin_mean12n(subj_id)
% M1-flat, regional linear obs, 12-node MEAN signal (no z-transform).
% Three fixes: detrend-only, default +1/128 A prior, per-node alpha.
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
TR = 0.5; n = 12;

d = load(fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', subj_id)), 'Y_mean');
signals = d.Y_mean;                         % T x 12
if size(signals,2) ~= n; signals = signals'; end
Y = detrend(signals);                       % preserve ΔF/F₀ baseline

A = ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
options.nonlinear=0; options.two_state=0; options.stochastic=0;
options.induced=1; options.centre=1; options.modality='Ca';
options.Fdcm=[0.01 0.2]; options.Nmax=32; options.maxit=256;
options.precision=log(16); options.verbose=0; options.linear_obs=0;
options.custom_g = @spm_gx_calcium_linear;

[pE,pC,x0] = spm_dcm_calcium_priors(A,B,C,D,options);
% Keep default +1/128 off-diagonal A prior (do NOT zero pE.A)

pE.alpha  = zeros(n,1);
pC.alpha  = ones(n,1) * 1/16;
pE.beta_y = zeros(n,1);
pC.beta_y = ones(n,1) * 1/64;
pE.Kd = 0; pC.Kd = 1/128;
pE.n  = 0; pC.n  = 1/128;

options.pE = pE; options.pC = pC;

DCM=struct(); DCM.a=A; DCM.b=B; DCM.c=C; DCM.d=D;
DCM.Y.y=Y; DCM.Y.dt=TR; DCM.Y.Q=[];
DCM.U.u=[]; DCM.U.dt=TR; DCM.M.nograph=1;
DCM.M.x=x0; DCM.options=options; DCM.options.maxnodes=max(16,n);
DCM.name=sprintf('DCM_sc_flat_reglin_mean12n_s%d', subj_id);

fprintf('S%d sc-flat RegLin mean12n... ', subj_id); rng(0); t0=tic;
DCM_est = spm_dcm_calcium_csd(DCM);
out = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', subj_id));
save(out,'DCM_est','-v7.3');
fprintf('done  F=%.2f  (%.1f min)\n', DCM_est.F, toc(t0)/60);
end
