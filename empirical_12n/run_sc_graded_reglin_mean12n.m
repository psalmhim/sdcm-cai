function run_sc_graded_reglin_mean12n(subj_id)
% M3 graded-SC prior, RegLin obs, 12-node Mean signal (k-means, no z-score).
% Variance for connected pairs = SC_norm(i,j)/16; unconnected = 1/128.
% Prior mean (pE.A) kept at spm_dcm_calcium_priors default (not zeroed).
%
% Saves: zebra/subject_N_DCM_sc_graded_reglin_mean12n.mat
% 2026 — IMAG-26-0111
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir  = './zebra/';
TR = 0.5; n = 12;
sc_idx_12 = [30, 66, 31, 23, 25, 17, 15, 67, 59, 61, 53, 51];
SC_data   = load('zebrafish_sc_mat.mat', 'sc_mat');
SC12      = double(SC_data.sc_mat(sc_idx_12, sc_idx_12));
for k = 1:n; SC12(k,k) = 0; end

d = load(fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', subj_id)), 'Y_mean');
signals = d.Y_mean;                         % T x 12
if size(signals,2) ~= n; signals = signals'; end
Y = detrend(signals);

A = ones(n,n);
B = zeros(n,n,0); C = zeros(n,0); D = zeros(n,n,0);
options.nonlinear=0; options.two_state=0; options.stochastic=0;
options.induced=1; options.centre=1; options.modality='Ca';
options.Fdcm=[0.01 0.2]; options.Nmax=32; options.maxit=256;
options.precision=log(16); options.verbose=0; options.linear_obs=0;
options.custom_g = @spm_gx_calcium_linear;

[pE,pC,x0] = spm_dcm_calcium_priors(A,B,C,D,options);

SC12_norm = SC12 / max(SC12(:));
for i = 1:n
    for j = 1:n
        if i == j; continue; end
        if SC12(i,j) > 0
            pC.A(i,j) = SC12_norm(i,j) * (1/16);
        else
            pC.A(i,j) = 1/128;
        end
    end
end

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
DCM.name=sprintf('DCM_sc_graded_reglin_mean12n_s%d', subj_id);

fprintf('S%d sc-graded RegLin-Mean 12n... ', subj_id); rng(0); t0=tic;
DCM_est = spm_dcm_calcium_csd(DCM);
out = fullfile(data_dir, sprintf('subject_%d_DCM_sc_graded_reglin_mean12n.mat', subj_id));
save(out,'DCM_est','-v7.3');
fprintf('done  F=%.2f  (%.1f min)\n', DCM_est.F, toc(t0)/60);
end
