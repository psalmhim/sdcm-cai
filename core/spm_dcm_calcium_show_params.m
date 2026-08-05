function spm_dcm_calcium_show_params(DCM, threshold)
%==========================================================================
% Display parameter estimates with labels for DCM calcium model
% FORMAT spm_dcm_calcium_show_params(DCM, threshold)
%
% DCM       - DCM structure with fields Ep, Cp, M.pE
% threshold - minimum abs(Ep-pE) to display (default: 0)
%
% Shows all parameters with their names, prior, posterior, and change
%==========================================================================
% Author: H.J. Park & ChatGPT, 2025
%==========================================================================

if nargin < 2, threshold = 0; end

pE = DCM.M.pE;
Ep = DCM.Ep;
Cp = DCM.Cp;

% Get parameter names and indices
[names, indices] = get_param_names(pE, DCM.n);

% Vectorize
vE = full(spm_vec(pE));
vEp = full(spm_vec(Ep));
vCp = full(diag(Cp));

% Calculate change from prior
delta = vEp - vE;
posterior_sd = sqrt(vCp);

fprintf('\n========================================================================\n');
fprintf('DCM Calcium Parameter Estimates\n');
fprintf('========================================================================\n');
fprintf('%-20s %6s %8s %10s %10s %10s\n', 'Parameter', 'Index', 'Prior', 'Posterior', 'Change', 'Post.SD');
fprintf('------------------------------------------------------------------------\n');

count_shown = 0;
count_total = length(delta);

for i = 1:length(delta)
    if abs(delta(i)) >= threshold
        fprintf('%s %6d %8.4f %10.4f %10.4f %10.4f\n', ...
            names{i}, i, vE(i), vEp(i), delta(i), posterior_sd(i));
        count_shown = count_shown + 1;
    end
end

fprintf('------------------------------------------------------------------------\n');
fprintf('Showing %d/%d parameters (threshold = %.4f)\n', count_shown, count_total, threshold);
fprintf('========================================================================\n\n');

% Plot all parameters with change from prior
figure('Name', 'All DCM Parameters', 'Position', [100 100 1200 600]);

subplot(2,2,1);
bar(delta);
xlabel('Parameter Index');
ylabel('Change from Prior');
title('All Parameters: Posterior - Prior');
set(gca,'xlim',[0 length(delta)+1]);
grid on;
hold on;
plot(xlim, [0 0], 'r--', 'LineWidth', 1);

subplot(2,2,2);
significant = abs(delta) > 0.01;
bar(find(significant), delta(significant));
xlabel('Parameter Index');
ylabel('Change from Prior');
title(sprintf('Significant Changes (|Δ| > 0.01): %d/%d params', sum(significant), length(delta)));
set(gca,'xlim',[0 length(delta)+1]);
grid on;
hold on;
plot(xlim, [0 0], 'r--', 'LineWidth', 1);

subplot(2,2,3);
errorbar(1:length(vEp), vEp, posterior_sd, 'o', 'MarkerSize', 3);
hold on;
plot(1:length(vE), vE, 'r.', 'MarkerSize', 8);
xlabel('Parameter Index');
ylabel('Value');
title('Prior (red) vs Posterior (blue) with SD');
legend('Posterior ± SD', 'Prior', 'Location', 'best');
grid on;

subplot(2,2,4);
% Group by parameter type
groups = get_param_groups(pE, DCM.n);
group_names = fieldnames(groups);
group_changes = zeros(length(group_names), 1);
for i = 1:length(group_names)
    idx = groups.(group_names{i});
    group_changes(i) = mean(abs(delta(idx)));
end
bar(group_changes);
set(gca, 'XTick', 1:length(group_names), 'XTickLabel', group_names, 'XTickLabelRotation', 45);
ylabel('Mean |Change|');
title('Average Change by Parameter Group');
grid on;

end

%==========================================================================
% Helper: Generate parameter names
%==========================================================================
function [names, indices] = get_param_names(pE, n)
names = {};
indices = struct();
idx = 1;

% A matrix
if isfield(pE, 'A')
    [na1, na2, na3] = size(pE.A);
    for k = 1:na3
        for i = 1:na1
            for j = 1:na2
                if na3 > 1
                    names{end+1} = sprintf('A(%d,%d,%d)', i, j, k);
                else
                    names{end+1} = sprintf('A(%d,%d)', i, j);
                end
                idx = idx + 1;
            end
        end
    end
end

% B matrix
if isfield(pE, 'B') && numel(pE.B) > 0
    [nb1, nb2, nb3] = size(pE.B);
    for k = 1:nb3
        for i = 1:nb1
            for j = 1:nb2
                names{end+1} = sprintf('B(%d,%d,%d)', i, j, k);
                idx = idx + 1;
            end
        end
    end
end

% C matrix
if isfield(pE, 'C') && numel(pE.C) > 0
    [nc1, nc2] = size(pE.C);
    for i = 1:nc1
        for j = 1:nc2
            names{end+1} = sprintf('C(%d,%d)', i, j);
            idx = idx + 1;
        end
    end
end

% D matrix
if isfield(pE, 'D') && numel(pE.D) > 0
    [nd1, nd2, nd3] = size(pE.D);
    for k = 1:nd3
        for i = 1:nd1
            for j = 1:nd2
                names{end+1} = sprintf('D(%d,%d,%d)', i, j, k);
                idx = idx + 1;
            end
        end
    end
end

% Calcium parameters (vector parameters)
vector_params = {'kappa_Ca', 'beta_Ca'};
for p = 1:length(vector_params)
    pname = vector_params{p};
    if isfield(pE, pname)
        nelem = numel(pE.(pname));
        for i = 1:nelem
            if nelem > 1
                names{end+1} = sprintf('%s(%d)', pname, i);
            else
                names{end+1} = pname;
            end
            idx = idx + 1;
        end
    end
end

% Scalar calcium parameters
scalar_params = {'R', 'tau_Ca', 'Kd', 'n', 'alpha', 'V0', 'beta_y'};
for p = 1:length(scalar_params)
    pname = scalar_params{p};
    if isfield(pE, pname)
        nelem = numel(pE.(pname));
        for i = 1:nelem
            if nelem > 1
                names{end+1} = sprintf('%s(%d)', pname, i);
            else
                names{end+1} = pname;
            end
            idx = idx + 1;
        end
    end
end

% Spectral parameters
if isfield(pE, 'a')
    na = numel(pE.a);
    for i = 1:na
        names{end+1} = sprintf('a(%d)', i);
        idx = idx + 1;
    end
end

if isfield(pE, 'b')
    nb = numel(pE.b);
    for i = 1:nb
        names{end+1} = sprintf('b(%d)', i);
        idx = idx + 1;
    end
end

if isfield(pE, 'c')
    nc = numel(pE.c);
    for i = 1:nc
        names{end+1} = sprintf('c(%d)', i);
        idx = idx + 1;
    end
end

indices = 1:length(names);
end

%==========================================================================
% Helper: Group parameters by type
%==========================================================================
function groups = get_param_groups(pE, n)
idx = 1;
groups = struct();

% Connectivity
if isfield(pE, 'A')
    nA = numel(pE.A);
    groups.Connectivity_A = idx:(idx+nA-1);
    idx = idx + nA;
end

if isfield(pE, 'B') && numel(pE.B) > 0
    nB = numel(pE.B);
    groups.Modulation_B = idx:(idx+nB-1);
    idx = idx + nB;
end

if isfield(pE, 'C') && numel(pE.C) > 0
    nC = numel(pE.C);
    groups.Input_C = idx:(idx+nC-1);
    idx = idx + nC;
end

if isfield(pE, 'D') && numel(pE.D) > 0
    nD = numel(pE.D);
    groups.Nonlinear_D = idx:(idx+nD-1);
    idx = idx + nD;
end

% Calcium dynamics
ca_params = {'kappa_Ca', 'beta_Ca', 'R', 'tau_Ca', 'Kd', 'n', 'alpha', 'V0', 'beta_y'};
for p = 1:length(ca_params)
    pname = ca_params{p};
    if isfield(pE, pname)
        nP = numel(pE.(pname));
        groups.(sprintf('Ca_%s', pname)) = idx:(idx+nP-1);
        idx = idx + nP;
    end
end

% Spectral
if isfield(pE, 'a')
    groups.Spectral_a = idx:(idx+numel(pE.a)-1);
    idx = idx + numel(pE.a);
end
if isfield(pE, 'b')
    groups.Spectral_b = idx:(idx+numel(pE.b)-1);
    idx = idx + numel(pE.b);
end
if isfield(pE, 'c')
    groups.Spectral_c = idx:(idx+numel(pE.c)-1);
    idx = idx + numel(pE.c);
end

end
