function param_map = spm_dcm_calcium_param_map(pE, n)
%==========================================================================
% Create a parameter index map for DCM calcium model
% FORMAT param_map = spm_dcm_calcium_param_map(pE, n)
%
% pE        - prior expectation structure
% n         - number of regions
%
% param_map - structure with parameter names and indices
%
% Returns a mapping of parameter names to their positions in spm_vec(pE)
% Useful for interpreting parameter plots and understanding which
% parameters are being estimated.
%==========================================================================
% Author: H.J. Park & ChatGPT, 2025
%==========================================================================

if nargin < 2
    if isfield(pE, 'A')
        n = size(pE.A, 1);
    else
        n = 1;
    end
end

param_map = struct();
param_map.names = {};
param_map.indices = [];
param_map.groups = {};
param_map.n_regions = n;

idx = 1;

% A matrix (connectivity)
if isfield(pE, 'A')
    [na1, na2, na3] = size(pE.A);
    start_idx = idx;
    for k = 1:na3
        for i = 1:na1
            for j = 1:na2
                if na3 > 1
                    param_map.names{end+1} = sprintf('A(%d,%d,%d)', i, j, k);
                else
                    param_map.names{end+1} = sprintf('A(%d,%d)', i, j);
                end
                param_map.indices(end+1) = idx;
                param_map.groups{end+1} = 'Connectivity';
                idx = idx + 1;
            end
        end
    end
    param_map.A_range = start_idx:(idx-1);
end

% B matrix (modulatory)
if isfield(pE, 'B') && numel(pE.B) > 0
    [nb1, nb2, nb3] = size(pE.B);
    start_idx = idx;
    for k = 1:nb3
        for i = 1:nb1
            for j = 1:nb2
                param_map.names{end+1} = sprintf('B(%d,%d,%d)', i, j, k);
                param_map.indices(end+1) = idx;
                param_map.groups{end+1} = 'Modulation';
                idx = idx + 1;
            end
        end
    end
    param_map.B_range = start_idx:(idx-1);
end

% C matrix (input)
if isfield(pE, 'C') && numel(pE.C) > 0
    [nc1, nc2] = size(pE.C);
    start_idx = idx;
    for i = 1:nc1
        for j = 1:nc2
            param_map.names{end+1} = sprintf('C(%d,%d)', i, j);
            param_map.indices(end+1) = idx;
            param_map.groups{end+1} = 'Input';
            idx = idx + 1;
        end
    end
    param_map.C_range = start_idx:(idx-1);
end

% D matrix (nonlinear)
if isfield(pE, 'D') && numel(pE.D) > 0
    [nd1, nd2, nd3] = size(pE.D);
    start_idx = idx;
    for k = 1:nd3
        for i = 1:nd1
            for j = 1:nd2
                param_map.names{end+1} = sprintf('D(%d,%d,%d)', i, j, k);
                param_map.indices(end+1) = idx;
                param_map.groups{end+1} = 'Nonlinear';
                idx = idx + 1;
            end
        end
    end
    param_map.D_range = start_idx:(idx-1);
end

% Calcium parameters - regional (per-region vectors)
regional_params = {'kappa_Ca', 'beta_Ca'};
for p = 1:length(regional_params)
    pname = regional_params{p};
    if isfield(pE, pname)
        nelem = numel(pE.(pname));
        start_idx = idx;
        for i = 1:nelem
            if nelem > 1
                param_map.names{end+1} = sprintf('%s(%d)', pname, i);
            else
                param_map.names{end+1} = pname;
            end
            param_map.indices(end+1) = idx;
            param_map.groups{end+1} = 'Calcium_regional';
            idx = idx + 1;
        end
        param_map.([pname '_range']) = start_idx:(idx-1);
    end
end

% Calcium parameters - global (scalar or shared)
global_params = {'R', 'tau_Ca', 'Kd', 'n', 'alpha', 'V0', 'beta_y'};
for p = 1:length(global_params)
    pname = global_params{p};
    if isfield(pE, pname)
        nelem = numel(pE.(pname));
        start_idx = idx;
        for i = 1:nelem
            if nelem > 1
                param_map.names{end+1} = sprintf('%s(%d)', pname, i);
            else
                param_map.names{end+1} = pname;
            end
            param_map.indices(end+1) = idx;
            param_map.groups{end+1} = 'Calcium_global';
            idx = idx + 1;
        end
        param_map.([pname '_range']) = start_idx:(idx-1);
    end
end

% Spectral parameters
spectral_params = {'a', 'b', 'c'};
for p = 1:length(spectral_params)
    pname = spectral_params{p};
    if isfield(pE, pname)
        nelem = numel(pE.(pname));
        start_idx = idx;
        for i = 1:nelem
            param_map.names{end+1} = sprintf('%s(%d)', pname, i);
            param_map.indices(end+1) = idx;
            param_map.groups{end+1} = 'Spectral';
            idx = idx + 1;
        end
        param_map.([pname '_range']) = start_idx:(idx-1);
    end
end

param_map.total_params = idx - 1;

% Create reverse lookup (index -> name)
param_map.index_to_name = containers.Map(param_map.indices, param_map.names);

% Summary by group
unique_groups = unique(param_map.groups);
param_map.summary = struct();
for i = 1:length(unique_groups)
    gname = unique_groups{i};
    param_map.summary.(gname) = find(strcmp(param_map.groups, gname));
end

end
