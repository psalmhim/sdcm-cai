function run_peb_iter2_collect_reglin_mean12n()
% Group PEB + BMR/BMA on RegLin-Mean iter2 12-node DCMs.
% Loads:  zebra/reglin_mean12n_iter2_subj_N.mat  (N=12:18)
% Saves:  zebra/PEB_iter2_sc_flat_reglin_mean12n_results.mat
%
% 2026 — IMAG-26-0111

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
n_nodes    = 12;
Np         = n_nodes * n_nodes;
node_names = {'lTeO_s','lTh','lP','lPT','lHb','lpRF', ...
              'rTeO_s','rTh','rP','rPT','rHb','rpRF'};

fprintf('=== RegLin-Mean iter2 group PEB (12-node) ===\n\n');

GCM_iter2 = {}; loaded_ids = [];
for s = subjects
    f = fullfile(data_dir, sprintf('reglin_mean12n_iter2_subj_%d.mat', s));
    if ~exist(f,'file'), fprintf('  Missing: S%d — skipping\n', s); continue; end
    d = load(f, 'DCM_est');
    GCM_iter2{end+1,1} = d.DCM_est;
    loaded_ids(end+1)  = s;
    fprintf('  S%d  F=%.2f\n', s, d.DCM_est.F);
end
Ns = numel(GCM_iter2);
fprintf('\n%d subjects loaded.\n\n', Ns);
if Ns < 3, error('Too few subjects (%d)', Ns); end

M.X      = ones(Ns,1);
M.Xnames = {'Group Mean'};
M.Q      = 'all';

fprintf('Running iter2 group PEB...\n');
PEB_iter2 = spm_dcm_peb(GCM_iter2, M, {'A'});
fprintf('PEB done. F=%.2f\n\n', PEB_iter2.F);

fprintf('Running BMR/BMA...\n');
p = gcp('nocreate'); if ~isempty(p), delete(p); end
[BMA_iter2, BMR_iter2] = spm_dcm_peb_bmc(PEB_iter2);
fprintf('BMA done.\n\n');

Ep_iter2 = reshape(full(real(BMA_iter2.Ep(1:Np,1))), n_nodes, n_nodes);
Pp_iter2 = reshape(full(real(BMA_iter2.Pp(1:Np))),   n_nodes, n_nodes);

Ep_peb  = reshape(full(real(PEB_iter2.Ep(1:Np,1))), n_nodes, n_nodes);
Cp_peb  = full(real(PEB_iter2.Cp(1:Np,1:Np)));
se_peb  = sqrt(abs(diag(Cp_peb)));
Pp_peb_mat = reshape(normcdf(abs(Ep_peb(:))./(se_peb+1e-12)), n_nodes, n_nodes);

fprintf('=== iter2 BMA (Pp>0.50) ===\n');
mask_od = ~eye(n_nodes,'logical');
[~,ord] = sort(Pp_iter2(:),'descend');
for idx = 1:numel(ord)
    li = ord(idx); [ri,ci] = ind2sub([n_nodes n_nodes],li);
    if ri==ci || Pp_iter2(ri,ci)<=0.50, break; end
    fprintf('  %s<-%s: Pp=%.3f  Ep=%+.4f\n', ...
        node_names{ri}, node_names{ci}, Pp_iter2(ri,ci), Ep_iter2(ri,ci));
end

for thr = [0.975 0.95 0.90]
    fprintf('Pp>%.3f: %d connections\n', thr, sum(Pp_iter2(mask_od)>thr));
end
fprintf('Pp_Gauss>0.975: %d connections\n', sum(Pp_peb_mat(mask_od)>0.975));

out_mat = fullfile(data_dir, 'PEB_iter2_sc_flat_reglin_mean12n_results.mat');
save(out_mat, 'PEB_iter2','BMA_iter2','BMR_iter2','GCM_iter2','loaded_ids', ...
    'node_names','Ep_iter2','Pp_iter2','Ep_peb','Pp_peb_mat','se_peb','-v7.3');
fprintf('\nSaved: %s\n', out_mat);
fprintf('=== iter2 PEB done ===\n');
end
