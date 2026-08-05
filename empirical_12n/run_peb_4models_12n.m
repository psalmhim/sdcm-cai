% PEB for all 4 obs models (M1-flat): RegHill-PC1, RegLin-PC1, RegHill-mean, RegLin-mean.
% Reports connections Pp>0.975 per model.
maxNumCompThreads(1);
home_dir = getenv('HOME');
addpath(fullfile(home_dir,'Dropbox/matlabwork/spm25'));
addpath(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
cd(fullfile(home_dir,'Dropbox/matlabwork/mnet0.92/dcmcai'));
spm('defaults','EEG'); spm_jobman('initcfg');

data_dir = './zebra/';
subjects  = 12:18;
n_subj    = numel(subjects);
n         = 12;

model_tags   = {'sc_flat_reghill_12n',    'sc_flat_reglin_12n', ...
                'sc_flat_reghill_mean12n', 'sc_flat_reglin_mean12n'};
model_labels = {'RegHill-PC1','RegLin-PC1','RegHill-mean','RegLin-mean'};

M_peb.Q = 'single';

results = struct();

for mi = 1:4
    tag   = model_tags{mi};
    label = model_labels{mi};
    fprintf('\n=== PEB: %s ===\n', label);

    GCMs = cell(n_subj,1);
    for si = 1:n_subj
        s  = subjects(si);
        fp = fullfile(data_dir, sprintf('subject_%d_DCM_%s.mat', s, tag));
        tmp = load(fp,'DCM_est');
        GCMs{si} = tmp.DCM_est;
    end

    [PEB, RCM] = spm_dcm_peb(GCMs, M_peb, {'A'});
    BMA = spm_dcm_peb_bmc(PEB);

    Ep_mat = full(reshape(BMA.Ep(1:n*n), n, n));
    Pp_mat = full(reshape(BMA.Pp(1:n*n), n, n));

    % Exclude diagonal (self-connections)
    off_diag = ~logical(eye(n));
    sig_mask = (Pp_mat > 0.975) & off_diag;
    n_sig = double(sum(sig_mask(:)));
    n_pos = double(sum(sig_mask(:) & (Ep_mat(:) > 0)));
    n_neg = double(sum(sig_mask(:) & (Ep_mat(:) < 0)));

    fprintf('  F_PEB = %.2f\n', PEB.F);
    fprintf('  Significant (Pp>0.975): %d total  (%d positive, %d negative)\n', n_sig, n_pos, n_neg);

    fprintf('  Connections (Pp>0.975):\n');
    [rows,cols] = find(sig_mask);
    for k = 1:numel(rows)
        i=rows(k); j=cols(k);
        fprintf('    A(%2d->%2d): Ep=%+.4f  Pp=%.3f\n', j, i, Ep_mat(i,j), Pp_mat(i,j));
    end

    results(mi).label  = label;
    results(mi).F_PEB  = PEB.F;
    results(mi).n_sig  = n_sig;
    results(mi).n_pos  = n_pos;
    results(mi).n_neg  = n_neg;
    results(mi).Ep_mat = Ep_mat;
    results(mi).Pp_mat = Pp_mat;
    results(mi).PEB    = PEB;
    results(mi).BMA    = BMA;

    save(fullfile(data_dir, sprintf('PEB_%s.mat', tag)), 'PEB','BMA','Ep_mat','Pp_mat','-v7.3');
    fprintf('  Saved PEB_%s.mat\n', tag);
end

fprintf('\n=== Summary ===\n');
fprintf('%-16s  %10s  %5s  %5s  %5s\n','Model','F_PEB','N_sig','N_pos','N_neg');
for mi = 1:4
    fprintf('%-16s  %10.2f  %5d  %5d  %5d\n', ...
        results(mi).label, results(mi).F_PEB, ...
        results(mi).n_sig, results(mi).n_pos, results(mi).n_neg);
end
fprintf('\nDone.\n');
