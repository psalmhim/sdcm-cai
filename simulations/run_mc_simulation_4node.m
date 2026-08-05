function Atrue = make_4node_groundtruth()

Atrue = [
   -0.5   0     -0.3   -0.1
    0.4  -0.5    0.2    0
    0     0.2   -0.5   -0.1
    0.1   0.3    0     -0.5
];

end


function sim = simulate_ca_dcm_spectral(cfg)

A  = cfg.A;
n  = size(A,1);

% ---------------------------
% Build minimal DCM structure
% ---------------------------
DCM = struct();
DCM.a = double(A~=0);
DCM.b = zeros(n,n,0);
DCM.c = zeros(n,0);
DCM.d = zeros(n,n,0);

% ====== 반드시 필요 ======
DCM.U.u    = zeros(n,1);      % no exogenous input
DCM.U.name = {'null'};
% ==========================

DCM.options.analysis  = 'CSD';
DCM.options.induced   = 1;
DCM.options.nonlinear = 0;
DCM.options.stochastic = 0;

% dummy Y (not used directly)
DCM.Y.y  = zeros(cfg.T,n);
DCM.Y.dt = cfg.TR;

% Priors
[pE,pC,x] = spm_dcm_calcium_priors( ...
    DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);

% overwrite with ground truth
pE.A        = A;
pE.tau_Ca  = log(cfg.tau);
pE.kappa_Ca   = log(cfg.kappa);
pE.R       = log(cfg.R);
pE.V0      = cfg.V0;
pE.beta_Ca = log(0.01*ones(n,1));

% Model specification
DCM.M.f  = @spm_fx_calcium;
DCM.M.g  = @spm_gx_calcium;
DCM.M.x  = x;
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.dt = cfg.TR;
DCM.M.N  = cfg.Nfreq;
DCM.M.Hz = linspace(0.02,1.0,cfg.Nfreq)';
DCM.M.u  = sparse(n,1);
DCM.M.p = 8;   % MAR model order (standard in spectral DCM)

% === THIS NOW WORKS ===
[Hs,Hz] = spm_csd_calcium_mtf(pE,DCM.M,DCM.U);

% add noise
noise = (1/cfg.SNR) * randn(size(Hs));

sim.Ycsd = Hs + noise;
sim.Hz   = Hz;
sim.Atrue = A;

end

function DCM = build_dcm_from_sim(sim)

n  = size(sim.Atrue,1);
nf = numel(sim.Hz);

DCM = struct();

% connectivity structure
DCM.a = double(sim.Atrue ~= 0);
DCM.b = zeros(n,n,0);
DCM.c = zeros(n,0);
DCM.d = zeros(n,n,0);

% ================================
% REQUIRED dummy time series
% ================================
DCM.Y.y  = zeros(512, n);   % dummy time series (never used)
DCM.Y.dt = 0.5;

% ================================
% spectral data
% ================================
DCM.Y.csd = sim.Ycsd;
DCM.Y.Hz  = sim.Hz;
DCM.Y.Q   = spm_dcm_csd_Q(sim.Ycsd);

% nuisance regressors (none)
DCM.Y.X0  = sparse(size(DCM.Y.Q,1),0);

% options
DCM.options.analysis = 'CSD';
DCM.options.induced  = 1;
DCM.options.order    = 8;
DCM.options.maxit    = 64;

end



function generate_sim_figures_4node(matfile)

load(matfile)

Atrue = Atrue;
n = size(Atrue,1);

Aest = cat(3,res.EpA);
Amean = mean(Aest,3);

% --- Figure 1: recovery ---
figure;
scatter(Atrue(~eye(n)),Amean(~eye(n)),'filled');
hold on;
plot([-0.5 0.5],[-0.5 0.5],'k--');
xlabel('A_{true}');
ylabel('A_{est}');
title('Connectivity recovery');

% --- Figure 2: ROC ---
truth = Atrue~=0 & ~eye(n);
pp = mean(cat(3,res.PpA),3);

[~,~,~,AUC] = perfcurve(truth(:),pp(:),1);
fprintf('AUC = %.3f\n',AUC);

figure;
perfcurve(truth(:),pp(:),1);
title(sprintf('Edge detection ROC (AUC=%.3f)',AUC));

end

addpath(fullfile(spm('Dir'),'toolbox','spectral'));
Nmc = 50;

Atrue = make_4node_groundtruth();

% ---- base parameters (read-only) ----
TR     = 0.5;
T      = 2400;
Nfreq  = 32;
tau    = 1.6;
kappa  = 5.5;
R      = 1.1;
V0     = -10;
SNR    = 10;

res(Nmc) = struct();   % preallocate

for i = 1:Nmc

    % ---- build cfg INSIDE parfor ----
    cfg = struct();
    cfg.A      = Atrue;
    cfg.TR     = TR;
    cfg.T      = T;
    cfg.Nfreq  = Nfreq;
    cfg.tau    = tau;
    cfg.kappa  = kappa;
    cfg.R      = R;
    cfg.V0     = V0;
    cfg.SNR    = SNR;

    % ---- simulate & invert ----
    sim = simulate_ca_dcm_spectral(cfg);
    DCM = build_dcm_from_sim(sim, cfg);
    DCM_est = spm_dcm_calcium_csd(DCM);

    % ---- store results ----
    res(i).EpA = DCM_est.Ep.A;
    res(i).PpA = DCM_est.Pp.A;
    res(i).F   = DCM_est.F;
end

save sim_results_4node.mat res Atrue

