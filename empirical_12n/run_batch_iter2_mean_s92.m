% Iter2 reinversion for RegLin-Mean subjects on server 92 (S15, S17).
% Requires: first-pass PEB at zebra/PEB_sc_flat_reglin_mean12n_results.mat
%           (generated on server 89 after all 7 flat RegLin-Mean DCMs done)
% 2026 — IMAG-26-0111
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath('/remotenas2/remotedata/share/matlabwork/spm25');
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

peb_file = './zebra/PEB_sc_flat_reglin_mean12n_results.mat';
fprintf('Waiting for first-pass PEB to appear via Dropbox sync...\n');
max_wait_min = 240;
for k = 1:max_wait_min
    if exist(peb_file,'file')
        d = load(peb_file);
        if isfield(d,'BMA')
            fprintf('PEB file found (%.0f min wait)\n', k);
            break;
        end
    end
    pause(60);
    if k == max_wait_min
        error('PEB file not found after %d min — run run_peb_sc_flat_reglin_mean12n first', max_wait_min);
    end
end

for s = [15 17]
    fprintf('\n=== Iter2 S%d ===\n', s);
    run_peb_iter2_single_reglin_mean12n(s);
end
fprintf('\n=== S92 iter2 DONE (S15, S17) ===\n');
