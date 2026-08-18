function run_peb_sc_flat_reglin_mean12n()
% First-pass group PEB + BMR/BMA on RegLin-Mean 12-node DCMs.
% Loads:  zebra/subject_N_DCM_sc_flat_reglin_mean12n.mat  (N=12:18)
% Saves:  zebra/PEB_sc_flat_reglin_mean12n_results.mat
%

maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath('/remotenas2/remotedata/share/matlabwork/spm25');
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG');
spm_jobman('initcfg');

data_dir   = './zebra/';
subjects   = 12:18;
node_names = {'lTeO_s','lTh','lP','lPT','lHb','lpRF', ...
              'rTeO_s','rTh','rP','rPT','rHb','rpRF'};
out_mat    = fullfile(data_dir, 'PEB_sc_flat_reglin_mean12n_results.mat');

fprintf('=== RegLin-Mean first-pass group PEB (12-node) ===\n');

GCM = {}; loaded_ids = [];
for s = subjects
    f = fullfile(data_dir, sprintf('subject_%d_DCM_sc_flat_reglin_mean12n.mat', s));
    if ~exist(f,'file'), fprintf('  Missing: S%d\n', s); continue; end
    d = load(f, 'DCM_est');
    GCM{end+1,1}     = d.DCM_est;
    loaded_ids(end+1) = s;
    fprintf('  S%d  F=%.2f\n', s, d.DCM_est.F);
end
Ns = numel(GCM);
fprintf('\n%d subjects loaded.\n\n', Ns);
if Ns < 2, error('Need >= 2 subjects'); end

M.X      = ones(Ns,1);
M.Xnames = {'Group Mean'};
M.Q      = 'all';

fprintf('Running PEB...\n');
PEB = spm_dcm_peb(GCM, M, {'A'});
fprintf('PEB done. F=%.2f\n\n', PEB.F);

save(out_mat, 'PEB', 'GCM', 'M', 'loaded_ids', 'node_names', '-v7.3');
fprintf('Interim saved: %s\n', out_mat);

fprintf('Running BMR/BMA...\n');
p = gcp('nocreate'); if ~isempty(p), delete(p); end
[BMA, BMR] = spm_dcm_peb_bmc(PEB);
fprintf('BMA done.\n\n');

save(out_mat, 'PEB', 'BMA', 'BMR', 'GCM', 'M', 'loaded_ids', 'node_names', '-v7.3');
fprintf('Saved: %s\n', out_mat);

n_nodes  = 12;
n_params = n_nodes^2;
Ep  = full(BMA.Ep(1:n_params,1));
Cp  = full(BMA.Cp(1:n_params,1:n_params));
pp  = 1 - spm_Ncdf(0, abs(Ep), sqrt(diag(Cp)));
fprintf('BMA: %d connections Pp>0.975\n', sum(pp>0.975));
fprintf('BMA: %d connections Pp>0.95\n',  sum(pp>0.95));
fprintf('\n=== First-pass PEB done ===\n');
end
