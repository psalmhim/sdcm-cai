function spm_dcm_calcium_results(DCM)
% Extended visualization for Calcium Spectral DCM
%
% FORMAT spm_dcm_calcium_results(DCM)
%
% Shows:
%   (1) Empirical vs predicted auto-spectra
%   (2) Empirical & predicted coherence matrices
%   (3) Pairwise connection-wise coherence spectra (frequency response)
%
% Author: HJ Park & ChatGPT, 2025
% -------------------------------------------------------------------------

if ischar(DCM)
    load(DCM, 'DCM');
end

if ~isfield(DCM, 'xY') || ~isfield(DCM.xY, 'g')
    error('DCM does not contain empirical spectra (xY.g)');
end

fprintf('\nPlotting Extended Calcium Spectral DCM Results...\n');

Hz = DCM.xY.Hz;
g_emp = DCM.xY.g;  % empirical cross-spectra
M  = DCM.M;
P  = DCM.Ep;

% Predict spectra using fitted parameters
[G_pred,~] = spm_csd_mtf_calcium(P,M,Hz);

n  = size(G_pred,1);
nf = length(Hz);

figure('Color','w','Name','Calcium Spectral DCM Results (Extended)', ...
       'Position',[100 100 1200 800]);

% -------------------------------------------------------------------------
% 1. Auto-spectra (Empirical vs Predicted)
% -------------------------------------------------------------------------
subplot(2,3,1)
hold on
for i = 1:n
    emp_auto = abs(squeeze(g_emp(:,i,i)));
    pred_auto = abs(squeeze(G_pred(i,i,:)));
    loglog(Hz, emp_auto, 'k-', 'LineWidth', 1.2)
    loglog(Hz, pred_auto, 'r--', 'LineWidth', 1.2)
end
xlabel('Frequency (Hz)')
ylabel('Power')
title('Auto-spectra (Empirical=Black, Predicted=Red)')
grid on
xlim([min(Hz) max(Hz)])

% -------------------------------------------------------------------------
% 2. Log Power Spectra (Linear frequency axis)
% -------------------------------------------------------------------------
subplot(2,3,2)
hold on
for i = 1:n
    emp_auto = abs(squeeze(g_emp(:,i,i)));
    pred_auto = abs(squeeze(G_pred(i,i,:)));
    plot(Hz, log10(emp_auto), 'k', 'LineWidth', 1.2)
    plot(Hz, log10(pred_auto), 'r--', 'LineWidth', 1.2)
end
xlabel('Frequency (Hz)')
ylabel('log_{10}(Power)')
title('Auto-spectra (log power)')
grid on

% -------------------------------------------------------------------------
% 3. Coherence Matrices
% -------------------------------------------------------------------------
subplot(2,3,3)
emp_coh = zeros(n,n);
pred_coh = zeros(n,n);

for i = 1:n
    for j = 1:n
        emp_coh(i,j) = mean(abs(squeeze(g_emp(:,i,j))) ./ ...
            sqrt(squeeze(g_emp(:,i,i)) .* squeeze(g_emp(:,j,j))));
        pred_coh(i,j) = mean(abs(squeeze(G_pred(i,j,:))) ./ ...
            sqrt(squeeze(G_pred(i,i,:)) .* squeeze(G_pred(j,j,:))));
    end
end

imagesc(emp_coh)
title('Empirical Coherence')
xlabel('Region')
ylabel('Region')
axis square
colorbar

subplot(2,3,6)
imagesc(pred_coh)
title('Predicted Coherence')
xlabel('Region')
ylabel('Region')
axis square
colorbar

% -------------------------------------------------------------------------
% 4. Pairwise Frequency Responses (Cross-Spectral Coherence)
% -------------------------------------------------------------------------
subplot(2,3,[4 5])
hold on

colors = lines(n);
legtxt = {};

for i = 1:n
    for j = i+1:n
        emp_spec = abs(squeeze(g_emp(:,i,j))) ./ ...
            sqrt(squeeze(g_emp(:,i,i)) .* squeeze(g_emp(:,j,j)));
        pred_spec = abs(squeeze(G_pred(i,j,:))) ./ ...
            sqrt(squeeze(G_pred(i,i,:)) .* squeeze(G_pred(j,j,:)));

        plot(Hz, emp_spec, '-', 'Color', colors(mod(i-1,size(colors,1))+1,:), ...
            'LineWidth', 1.2, 'DisplayName', sprintf('Emp %d-%d',i,j));
        plot(Hz, pred_spec, '--', 'Color', colors(mod(i-1,size(colors,1))+1,:), ...
            'LineWidth', 1.2, 'HandleVisibility','off');
        legtxt{end+1} = sprintf('%d-%d', i, j);
    end
end

xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Pairwise Coherence Spectra (Emp=Solid, Pred=Dashed)')
legend(legtxt, 'Location', 'northeastoutside')
grid on
xlim([min(Hz) max(Hz)])

sgtitle('Spectral DCM for Calcium Imaging: Empirical vs Predicted','FontWeight','bold')

fprintf('Done.\n\n');
fprintf('Model Free Energy: %.3f\n', DCM.F);
disp('Posterior connectivity (A-matrix):')
disp(DCM.Ep.a);
end
