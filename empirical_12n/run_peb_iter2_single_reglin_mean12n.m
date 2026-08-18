function run_peb_iter2_single_reglin_mean12n(subj_id)
% Iter2 re-inversion for RegLin-Mean 12-node model.
% Uses empirical Bayes prior from first-pass group PEB.
%
% Loads:  zebra/PEB_sc_flat_reglin_mean12n_results.mat  (group prior)
%         zebra/subject_N_DCM_sc_flat_reglin_mean12n.mat (original DCM)
%         zebra/subject_N_meantrace_12n.mat              (Mean signals)
% Saves:  zebra/reglin_mean12n_iter2_subj_N.mat
%

maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath('/remotenas2/remotedata/share/matlabwork/spm25');
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG');
spm_jobman('initcfg');

data_dir = './zebra/';
n_nodes  = 12;
Np       = n_nodes * n_nodes;

%% 1) Load first-pass group PEB
peb_file = fullfile(data_dir, 'PEB_sc_flat_reglin_mean12n_results.mat');
if ~exist(peb_file,'file'), error('Run run_peb_sc_flat_reglin_mean12n first.'); end
peb_data = load(peb_file, 'PEB', 'GCM');
PEB = peb_data.PEB;
GCM = peb_data.GCM;

%% 2) Extract empirical Bayes prior for A
A_group = reshape(full(real(PEB.Ep(1:Np, 1))), n_nodes, n_nodes);
fprintf('S%d — group mean |A_offdiag|=%.4f\n', subj_id, ...
    mean(abs(A_group(~logical(eye(n_nodes))))));

pC_A_vec = [];
if isfield(PEB,'Ce') && ~isempty(PEB.Ce)
    Ce_diag = diag(full(real(PEB.Ce)));
    if numel(Ce_diag) >= Np
        pC_A_vec = Ce_diag(1:Np);
    end
end
if isempty(pC_A_vec)
    A_stack = zeros(n_nodes, n_nodes, numel(GCM));
    for i = 1:numel(GCM)
        A_stack(:,:,i) = GCM{i}.Ep.A;
    end
    pC_A_vec = reshape(var(A_stack, 0, 3), Np, 1);
end
pC_A_vec = max(pC_A_vec, 1e-4);

%% 3) Load original DCM to inherit exact options/priors
orig_file = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', subj_id));
orig      = load(orig_file, 'DCM_est');
DCM_orig  = orig.DCM_est;

%% 4) Load Mean signals
sig_file = fullfile(data_dir, sprintf('subject_%d_meantrace_12n.mat', subj_id));
d = load(sig_file, 'Y_mean');
signals = d.Y_mean;                         % T x 12
if size(signals,2) ~= n_nodes, signals = signals'; end
Y = detrend(signals);

%% 5) Build iter2 DCM — inherit options, override only pE.A / pC.A
pE = DCM_orig.M.pE;
pC = DCM_orig.M.pC;
pE.A = A_group;
pC.A = reshape(pC_A_vec, n_nodes, n_nodes);

DCM           = struct();
DCM.a         = DCM_orig.a;
DCM.b         = DCM_orig.b;
DCM.c         = DCM_orig.c;
DCM.d         = DCM_orig.d;
DCM.Y.y       = Y;
DCM.Y.dt      = DCM_orig.Y.dt;
DCM.Y.Q       = [];
DCM.U.u       = [];
DCM.U.dt      = DCM_orig.U.dt;
DCM.M         = DCM_orig.M;
DCM.M.pE      = pE;
DCM.M.pC      = pC;
DCM.options   = DCM_orig.options;
DCM.options.pE = pE;
DCM.options.pC = pC;
DCM.options.maxnodes = max(16, n_nodes);
DCM.name      = sprintf('reglin_mean12n_iter2_subj_%d', subj_id);

fprintf('Inverting S%d RegLin-Mean iter2... ', subj_id);
rng(0); t0 = tic;
try
    DCM_est = spm_dcm_calcium_csd(DCM);
    out = fullfile(data_dir, sprintf('reglin_mean12n_iter2_subj_%d.mat', subj_id));
    save(out, 'DCM_est', 'A_group', 'pC_A_vec', '-v7.3');
    fprintf('done  F=%.2f  (%.1f min)\n', DCM_est.F, toc(t0)/60);
catch ME
    fprintf('FAILED: %s\n', ME.message);
end
end
