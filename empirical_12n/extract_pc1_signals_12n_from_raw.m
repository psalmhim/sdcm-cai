function extract_pc1_signals_12n_from_raw(subj_id, w_sp)
% Extract 12-node PC1-active signals directly from raw calcium data.
%
% Node ordering (6 bilateral pairs):
%   lTeO_s, lTh, lP, lPT, lHb, lpRF
%   rTeO_s, rTh, rP, rPT, rHb, rpRF
%
% Pipeline per node:
%   1. Load cells from atlas ROI
%   2. Select active cells (variance > median)
%   For TeO_s only:
%     2a. Joint spatial+temporal k-means (k=2) on active cells
%         using teo_kmeans_split(X_active, xyz_n_active, w_sp)
%         w_sp: weight for spatial features [0=pure temporal, 1=pure spatial]
%     2b. Superficial cluster identified by lower mean DV
%   3. Cell-level standardisation:
%      detrend → mean-centre → unit-variance each cell
%   4. Population signal = PC1 of standardised active cells (PC1-active)
%   5. Sign-align to ROI mean
%
% Subject-level global scaling is applied in DCM scripts, NOT here.
%
% INPUTS
%   subj_id  subject index (12-18)
%   w_sp     spatial weight for TeO k-means split (default 0.5)
%            0 = pure time-series clustering
%            1 = pure spatial PCA split (legacy behaviour)
%
% OUTPUT: zebra/subject_N_pc1signals_12n.mat
%   pc1_active_signals  [12 x T]
%   node_names, node_idx
%
% 2026 — IMAG-26-0111

if nargin < 2 || isempty(w_sp), w_sp = 0.5; end

maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));

raw_dir  = '/remotenas2/remotedata/share/zebrafish/original_raw_data';
data_dir = fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai/zebra');
min_cells = 10;

fname = fullfile(raw_dir, sprintf('subject_%d_data.mat', subj_id));
if ~exist(fname,'file')
    error('Raw data not found: %s', fname);
end

fprintf('Subject %d (12-node direct extraction)...\n', subj_id);
d          = load(fname, 's_data', 'result_idx', 'CellXYZ');
s_data     = d.s_data;
result_idx = d.result_idx;
cell_xyz   = d.CellXYZ;   % [N_cells x 3]
T          = size(s_data, 2);

%% Node definitions
node_names = {'lTeO_s','lTh','lP','lPT','lHb','lpRF', ...
              'rTeO_s','rTh','rP','rPT','rHb','rpRF'};
node_idx   = [30 31 23 25 17 15  66 67 59 61 53 51];
n_nodes    = 12;

pc1_active_signals = nan(n_nodes, T);

%% --- Non-TeO nodes (slots 2-6 left, 8-12 right) ---
nonTeO_atlas = [31 23 25 17 15   67 59 61 53 51];
nonTeO_slots = [2  3  4  5  6    8  9  10 11 12];

for k = 1:numel(nonTeO_atlas)
    ri  = nonTeO_atlas(k);
    idx = result_idx{ri};
    if numel(idx) < min_cells; continue; end

    X = s_data(idx, :);                          % [n_cells x T]
    X(any(~isfinite(X),2), :) = [];
    if size(X,1) < min_cells; continue; end

    roi_mean = mean(X, 1);                        % [1 x T] for sign-align

    % Active-cell selection (top 50% by variance)
    cell_var = var(X, 0, 2);
    Xa = X(cell_var > median(cell_var), :);
    if size(Xa,1) < min_cells; continue; end

    % Cell-level standardisation → mean of standardised cells
    Xa_n = detrend(Xa')';                         % detrend each cell [n_cells x T]
    Xa_n = Xa_n - mean(Xa_n, 2);                 % mean-centre each cell
    cs   = std(Xa_n, 0, 2); cs(cs < 1e-10) = 1e-10;
    Xa_n = Xa_n ./ cs;                            % unit-variance each cell

    [~, sc_] = pca(Xa_n', 'Centered', false);
    pop_sig  = sc_(:,1)';                         % PC1 scores [1 x T]
    if corr(roi_mean', pop_sig') < 0; pop_sig = -pop_sig; end

    pc1_active_signals(nonTeO_slots(k), :) = pop_sig;
    fprintf('  %s: %d active cells -> PC1\n', node_names{nonTeO_slots(k)}, size(Xa,1));
end

%% --- TeO superficial nodes (slots 1 and 7) ---
teo_defs = struct('atlas_ri', {30, 66}, 'slot', {1, 7}, ...
                  'side_name', {'L', 'R'}, 'mirror_lr', {false, true});

for td = 1:2
    ri        = teo_defs(td).atlas_ri;
    slot_s    = teo_defs(td).slot;
    side_name = teo_defs(td).side_name;
    mirror_lr = teo_defs(td).mirror_lr;

    idx = result_idx{ri};
    if numel(idx) < 2*min_cells
        fprintf('  %sTeO: too few cells (%d), skipping\n', side_name, numel(idx));
        continue;
    end

    X_all    = s_data(idx, :);
    fin_mask = all(isfinite(X_all), 2);      % identify non-finite rows first
    X        = X_all(fin_mask, :);           % filter X and xyz together
    raw_xyz  = cell_xyz(idx, :);
    raw_xyz  = raw_xyz(fin_mask, :);         % same mask, same row count as X

    nc = size(X, 1);
    if nc < 2*min_cells; continue; end

    roi_mean = mean(X, 1);

    %% Step 1: Normalise cell xyz to [0,1] bounding box
    lo_xyz  = min(raw_xyz);
    hi_xyz  = max(raw_xyz);
    rng_xyz = max(hi_xyz - lo_xyz, 1e-9);
    xyz_n   = (raw_xyz - lo_xyz) ./ rng_xyz;

    if mirror_lr
        xyz_n(:, 2) = 1 - xyz_n(:, 2);
    end

    %% Step 2: Active-cell mask (top 50% variance over ALL cells)
    cell_var = var(X, 0, 2);
    act_mask = cell_var > median(cell_var);
    n_act    = sum(act_mask);

    fprintf('  %sTeO: %d total, %d active\n', side_name, nc, n_act);

    if n_act < 2*min_cells
        fprintf('  %sTeO: too few active cells (%d), skipping\n', side_name, n_act);
        continue;
    end

    Xa     = X(act_mask, :);           % [n_act x T]
    xyz_na = xyz_n(act_mask, :);       % [n_act x 3]

    %% Step 3: Joint spatial+temporal k-means split (k=2) on active cells
    km_opts.n_ts_pc  = 15;
    km_opts.n_rep    = 50;
    km_opts.rng_seed = 42 + subj_id;  % reproducible per subject
    km_opts.verbose  = true;

    [lbl_deep_act, km_info] = teo_kmeans_split(Xa, xyz_na, w_sp, km_opts);
    fprintf('  %sTeO: w_sp=%.2f  superf=%d(DV=%.3f,r=%.3f)  deep=%d(DV=%.3f,r=%.3f)\n', ...
            side_name, w_sp, ...
            km_info.n_superf, km_info.dv_superf, km_info.r_superf, ...
            km_info.n_deep,   km_info.dv_deep,   km_info.r_deep);

    mask_superf_active = ~lbl_deep_act;   % [n_act x 1] logical
    n_s = sum(mask_superf_active);

    if n_s < min_cells
        fprintf('  %sTeO_s: too few cells in superficial cluster (%d), skipping\n', side_name, n_s);
        continue;
    end

    Xs = Xa(mask_superf_active, :);       % [n_s x T]

    %% Step 4: Cell-level standardisation → PC1 of superficial cluster
    Xs_n = detrend(Xs')';                 % detrend each cell
    Xs_n = Xs_n - mean(Xs_n, 2);         % mean-centre
    cs_s = std(Xs_n, 0, 2); cs_s(cs_s < 1e-10) = 1e-10;
    Xs_n = Xs_n ./ cs_s;                  % unit-variance

    %% Step 5: Population signal = PC1
    [~, sc_s] = pca(Xs_n', 'Centered', false);
    pc1s      = sc_s(:,1)';               % PC1 scores [1 x T]
    if corr(roi_mean', pc1s') < 0; pc1s = -pc1s; end

    pc1_active_signals(slot_s, :) = pc1s;
    fprintf('    %sTeO_s: %d cells -> PC1 extracted\n', side_name, n_s);
end

%% Save
out_file = fullfile(data_dir, sprintf('subject_%d_pc1signals_12n.mat', subj_id));
save(out_file, 'pc1_active_signals', 'node_names', 'node_idx', '-v7.3');
fprintf('Saved: %s  [12 x %d]\n', out_file, T);

% Quick check
r = corrcoef(pc1_active_signals');
mask_od = ~eye(12,'logical');
fprintf('  mean|r|=%.3f  lTeOs<->rTeOs=%.3f\n', ...
    mean(abs(r(mask_od))), r(1,7));
end
