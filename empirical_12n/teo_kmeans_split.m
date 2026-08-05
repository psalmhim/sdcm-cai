function [lbl_deep, info] = teo_kmeans_split(X, xyz_n, w_sp, opts)
%==========================================================================
% Split TeO cells into superficial / deep subpopulations using joint
% spatial + temporal k-means (k=2).
%
% FORMAT [lbl_deep, info] = teo_kmeans_split(X, xyz_n, w_sp, opts)
%
% INPUTS
%   X      [n_cells x T]  cell time series (raw or z-scored)
%   xyz_n  [n_cells x 3]  normalised cell xyz, each axis in [0,1]
%   w_sp   scalar in [0,1] weight for spatial features (1=pure spatial,
%                          0=pure temporal).  Default 0.5.
%   opts   struct (optional):
%     .n_ts_pc   number of time series PCs to retain (default 15)
%     .n_rep     k-means replicates (default 50)
%     .rng_seed  RNG seed (default 42)
%     .verbose   print diagnostic info (default true)
%
% OUTPUT
%   lbl_deep  [n_cells x 1] logical: true = deep, false = superficial
%             Convention: superficial cluster has lower mean DV (xyz col 3).
%   info      struct with cluster diagnostics
%
% DISTANCE METRIC
%   Features are split into two blocks:
%     spatial : xyz_n (3 columns), std-normalised across cells
%     temporal: first n_ts_pc PCs of z-scored time series, std-normalised
%
%   Weighted concatenation:
%     feat = [sqrt(w_sp)*(spatial_block), sqrt(1-w_sp)*(temporal_block)]
%
%   K-means with Euclidean distance on this combined feature matrix.
%   Euclidean distance on z-scored time series is proportional to
%   sqrt(2*(1-corr)), so the temporal block implements a correlation-like
%   similarity weighted by spatial proximity.
%
% 2026 - IMAG-26-0111 (joint spatial-temporal TeO split)
%==========================================================================

if nargin < 3 || isempty(w_sp),  w_sp = 0.5;  end
if nargin < 4,                    opts = struct(); end

try, n_ts_pc  = opts.n_ts_pc;  catch, n_ts_pc  = 15;   end
try, n_rep    = opts.n_rep;    catch, n_rep    = 50;   end
try, seed     = opts.rng_seed; catch, seed     = 42;   end
try, verbose  = opts.verbose;  catch, verbose  = true; end

[n_cells, T] = size(X);

%--------------------------------------------------------------------------
% (1) Temporal features: z-score each cell, PCA-reduce
%--------------------------------------------------------------------------
Xz = detrend(X')';                    % detrend per cell [n_cells x T]
Xz = Xz - mean(Xz, 2);              % mean-centre
cs = std(Xz, 0, 2); cs(cs < 1e-10) = 1e-10;
Xz = Xz ./ cs;                        % unit-variance per cell

n_comp = min(n_ts_pc, min(n_cells, T) - 1);
[~, ts_scores] = pca(Xz, 'NumComponents', n_comp, 'Centered', false);
% ts_scores: [n_cells x n_comp]

%--------------------------------------------------------------------------
% (2) Spatial features: use xyz_n as-is (already in [0,1])
%--------------------------------------------------------------------------
sp_feat = xyz_n;    % [n_cells x 3]

%--------------------------------------------------------------------------
% (3) Std-normalise each feature block
%--------------------------------------------------------------------------
sp_std  = std(sp_feat, 0, 1);
sp_std(sp_std < 1e-10) = 1e-10;
sp_norm = sp_feat ./ sp_std;           % [n_cells x 3], unit column-std

ts_std  = std(ts_scores, 0, 1);
ts_std(ts_std < 1e-10) = 1e-10;
ts_norm = ts_scores ./ ts_std;         % [n_cells x n_comp], unit column-std

%--------------------------------------------------------------------------
% (4) Weighted concatenation
%--------------------------------------------------------------------------
feat = [sqrt(w_sp) * sp_norm, sqrt(1 - w_sp) * ts_norm];

%--------------------------------------------------------------------------
% (5) K-means (k=2)
%--------------------------------------------------------------------------
rng(seed);
[lbl, C, sumd] = kmeans(feat, 2, ...
    'Distance',    'sqeuclidean', ...
    'Replicates',  n_rep, ...
    'MaxIter',     500, ...
    'Display',     'off');
% lbl: [n_cells x 1], values 1 or 2

%--------------------------------------------------------------------------
% (6) Assign superficial / deep by DV convention
%     Superficial = lower mean DV (xyz column 3)
%--------------------------------------------------------------------------
dv_col = 3;
dv_c1  = mean(xyz_n(lbl == 1, dv_col));
dv_c2  = mean(xyz_n(lbl == 2, dv_col));

if dv_c1 < dv_c2
    % cluster 1 is superficial (lower DV)
    lbl_deep = lbl == 2;
else
    % cluster 2 is superficial (lower DV)
    lbl_deep = lbl == 1;
end

%--------------------------------------------------------------------------
% (7) Diagnostics
%--------------------------------------------------------------------------
n_s = sum(~lbl_deep);
n_d = sum( lbl_deep);

% Mean within-cluster correlation for temporal block
r_s = mean_pairwise_corr(X(~lbl_deep, :));
r_d = mean_pairwise_corr(X( lbl_deep, :));

% Spatial separation (DV)
dv_s = mean(xyz_n(~lbl_deep, dv_col));
dv_d = mean(xyz_n( lbl_deep, dv_col));

% Inertia ratio (between / total)
total_inertia = sum(sumd);

info = struct('n_superf', n_s, 'n_deep', n_d, ...
              'dv_superf', dv_s, 'dv_deep', dv_d, ...
              'r_superf', r_s, 'r_deep', r_d, ...
              'kmeans_inertia', total_inertia, ...
              'w_sp', w_sp, 'n_ts_pc', n_comp, ...
              'centroid_sp', C(:, 1:3), ...
              'centroid_ts', C(:, 4:end));

if verbose
    fprintf('    K-means split: superf=%d(DV=%.3f, r=%.3f)  deep=%d(DV=%.3f, r=%.3f)\n', ...
            n_s, dv_s, r_s, n_d, dv_d, r_d);
    fprintf('    w_sp=%.2f  n_ts_pc=%d  inertia=%.2f\n', w_sp, n_comp, total_inertia);
end

end

%--------------------------------------------------------------------------
function r = mean_pairwise_corr(X)
% Mean pairwise correlation for a sample of cells [fast approximation]
if size(X, 1) < 2
    r = NaN; return
end
n = min(size(X, 1), 200);   % subsample if large for speed
if size(X, 1) > n
    idx = randperm(size(X, 1), n);
    X = X(idx, :);
end
Xz = X - mean(X, 2);
nrm = sqrt(sum(Xz.^2, 2));
nrm(nrm < 1e-12) = 1e-12;
Xz = Xz ./ nrm;
R  = Xz * Xz';
mask = ~eye(size(R,1), 'logical');
r  = mean(R(mask));
end
