function spm_dcm_calcium_results_summary(DCM)
% Extended visualization and quantitative summary for Calcium Spectral DCM
%
% FORMAT spm_dcm_calcium_results_summary(DCM)
%
% Plots:
%   (1) Empirical vs predicted auto-spectra
%   (2) Empirical & predicted coherence matrices
%   (3) Pairwise coherence spectra (empirical/predicted)
%
% Prints:
%   (4) Summary table of mean coherence and power in low/mid/high bands
%
% Author: HJ Park & ChatGPT, 2025
% -------------------------------------------------------------------------

if ischar(DCM)
    load(DCM,'DCM');
end
if ~isfield(DCM,'xY') 
    DCM.xY.g  = DCM.Y.csd;       % empirical cross spectra
    DCM.xY.Hz = DCM.Y.Hz;        % frequency vector
    DCM.xY.name = arrayfun(@(i) sprintf('Region %d', i), 1:DCM.n, 'UniformOutput', false);
end

Hz = DCM.Y.Hz(:);
g_emp = DCM.Y.csd;
M  = DCM.M;
P  = DCM.Ep;

% Predicted spectra
[G_pred,~] = spm_csd_calcium_mtf(P,M,Hz);

n  = size(G_pred,3);
fprintf('\n=== Spectral DCM for Calcium Imaging: Results Summary ===\n');
fprintf('Model Free Energy: %.3f\n\n', DCM.F);

% -------------------------------------------------------------------------
% PLOTTING
% -------------------------------------------------------------------------
figure('Color','w','Name','Calcium Spectral DCM Results (Summary)',...
       'Position',[100 100 1200 800]);

% (1) Auto-spectra (log-log)
subplot(2,3,1); hold on
for i = 1:n
    emp_auto = abs(squeeze(g_emp(:,i,i)));
    pred_auto = abs(squeeze(G_pred(:,i,i)));
    loglog(Hz, emp_auto, 'k-', 'LineWidth', 1.3);
    loglog(Hz, pred_auto, 'r--', 'LineWidth', 1.3);
end
xlabel('Frequency (Hz)'); ylabel('Power');
title('Auto-spectra (Emp=Black, Pred=Red)'); grid on

% (2) Log power spectra
subplot(2,3,2); hold on
for i = 1:n
    emp_auto = abs(squeeze(g_emp(:,i,i)));
    pred_auto = abs(squeeze(G_pred(:,i,i)));
    plot(Hz, log10(emp_auto), 'k', 'LineWidth', 1.3);
    plot(Hz, log10(pred_auto), 'r--', 'LineWidth', 1.3);
end
xlabel('Frequency (Hz)'); ylabel('log_{10}(Power)');
title('Auto-spectra (linear freq)'); grid on

% (3,6) Coherence matrices
subplot(2,3,3)
emp_coh = zeros(n,n); pred_coh = zeros(n,n);
for i = 1:n
    for j = 1:n
        emp_coh(i,j) = mean(abs(squeeze(g_emp(:,i,j))) ./ ...
            sqrt(squeeze(g_emp(:,i,i)).*squeeze(g_emp(:,j,j))));
        pred_coh(i,j)= mean(abs(squeeze(G_pred(:,i,j))) ./ ...
            sqrt(squeeze(G_pred(:,i,i)).*squeeze(G_pred(:,j,j))));
    end
end
imagesc(abs(emp_coh)); title('abs Empirical Coherence');
xlabel('Region'); ylabel('Region'); axis square; colorbar

subplot(2,3,6)
imagesc(abs(pred_coh)); title('abs Predicted Coherence');
xlabel('Region'); ylabel('Region'); axis square; colorbar

% (4-5) Pairwise coherence spectra
subplot(2,3,[4 5]); hold on
colors = lines(n);
for i = 1:n
    for j = i+1:n
        emp_spec = abs(squeeze(g_emp(:,i,j))) ./ ...
            sqrt(squeeze(g_emp(:,i,i)) .* squeeze(g_emp(:,j,j)));
        pred_spec = abs(squeeze(G_pred(:,i,j))) ./ ...
            sqrt(squeeze(G_pred(i,i)) .* squeeze(G_pred(:,j,j)));
        plot(Hz, emp_spec, '-', 'Color', colors(mod(i-1,size(colors,1))+1,:), 'LineWidth', 1.3);
        plot(Hz, pred_spec, '--', 'Color', colors(mod(i-1,size(colors,1))+1,:), 'LineWidth', 1.3);
    end
end
xlabel('Frequency (Hz)'); ylabel('Coherence');
title('Pairwise Coherence (Emp=Solid, Pred=Dashed)');
grid on; xlim([min(Hz) max(Hz)]);
sgtitle('Spectral DCM for Calcium Imaging: Model Fit','FontWeight','bold');

% -------------------------------------------------------------------------
% SUMMARY STATISTICS BY FREQUENCY BAND
% -------------------------------------------------------------------------
bands = [0.1 0.5; 0.5 2; 2 5];   % Hz ranges: low / mid / high
bandnames = {'Low','Mid','High'};
fprintf('Frequency Bands (Hz):\n');
for b = 1:size(bands,1)
    fprintf('  %s: %.1f - %.1f Hz\n', bandnames{b}, bands(b,1), bands(b,2));
end

% mean power per region, mean coherence per connection
fprintf('\nMean Power & Coherence per Band:\n');
fprintf('------------------------------------------------------------\n');
fprintf('%10s %8s %8s %8s %8s %8s\n','Region/Conn','Type','Low','Mid','High','Overall');
fprintf('------------------------------------------------------------\n');

% Auto-spectral power per ROI
for i = 1:n
    P_auto = abs(squeeze(G_pred(:,i,i)));
    for b = 1:size(bands,1)
        idx = Hz >= bands(b,1) & Hz < bands(b,2);
        bandpow(b) = mean(P_auto(idx));
    end
    totalpow = mean(P_auto);
    fprintf('%10s %8s %8.3f %8.3f %8.3f %8.3f\n',...
        sprintf('ROI%d',i),'Power',bandpow(1),bandpow(2),bandpow(3),totalpow);
end

% Coherence per connection
for i = 1:n
    for j = i+1:n
        Cxy = abs(squeeze(G_pred(:,i,j))) ./ ...
            sqrt(squeeze(G_pred(:,i,i)).*squeeze(G_pred(:,j,j)));
        for b = 1:size(bands,1)
            idx = Hz >= bands(b,1) & Hz < bands(b,2);
            bandcoh(b) = mean(Cxy(idx));
        end
        totalcoh = mean(Cxy);
        fprintf('%10s %8s %8.3f %8.3f %8.3f %8.3f\n',...
            sprintf('%d-%d',i,j),'Coher',bandcoh(1),bandcoh(2),bandcoh(3),totalcoh);
    end
end
fprintf('------------------------------------------------------------\n');
fprintf('All values are averages of predicted spectra over each band.\n');
fprintf('------------------------------------------------------------\n\n');

end
