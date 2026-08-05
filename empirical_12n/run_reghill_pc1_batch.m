% Batch re-run of RegHill-PC1 (12n) with bug-fixed spm_dcm_mtf_calcium.m
% Bug fix (Jul 15): operating point no longer forced to Ca=0.5 when x(:,1)=0
% This overwrites the buggy Jul-14 DCMs.
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

subjects = [12 13 14 15 16 17 18];
data_dir = './zebra/';
TR = 0.5; n = 12;

A = ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
options.nonlinear=0; options.two_state=0; options.stochastic=0;
options.induced=1; options.centre=1; options.modality='Ca';
options.Fdcm=[0.01 0.2]; options.Nmax=32; options.maxit=256;
options.precision=log(16); options.verbose=0; options.linear_obs=0;

[pE,pC,x0] = spm_dcm_calcium_priors(A,B,C,D,options);

pE.alpha  = zeros(n,1);
pC.alpha  = ones(n,1) * 1/16;
pE.beta_y = zeros(n,1);
pC.beta_y = ones(n,1) * 1/64;
pE.Kd = 0; pC.Kd = 1/128;
pE.n  = 0; pC.n  = 1/128;

options.pE = pE; options.pC = pC;

for si = 1:numel(subjects)
    s = subjects(si);
    fprintf('\n=== Subject %d (%d/%d) ===\n', s, si, numel(subjects));

    d = load(fullfile(data_dir, sprintf('subject_%d_pc1signals_12n.mat', s)), 'pc1_active_signals');
    signals = d.pc1_active_signals;
    if size(signals,1) > size(signals,2); signals = signals'; end
    Y = detrend(signals');

    DCM = struct();
    DCM.a=A; DCM.b=B; DCM.c=C; DCM.d=D;
    DCM.Y.y=Y; DCM.Y.dt=TR; DCM.Y.Q=[];
    DCM.U.u=[]; DCM.U.dt=TR; DCM.M.nograph=1;
    DCM.M.x=x0; DCM.options=options; DCM.options.maxnodes=max(16,n);
    DCM.name=sprintf('DCM_sc_flat_reghill_12n_s%d', s);

    rng(0); t0 = tic;
    DCM_est = spm_dcm_calcium_csd(DCM);
    out = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reghill_12n.mat', s));
    save(out,'DCM_est','-v7.3');
    fprintf('S%d done  F=%.2f  (%.1f min)\n', s, DCM_est.F, toc(t0)/60);
end

fprintf('\n=== All subjects done ===\n');
