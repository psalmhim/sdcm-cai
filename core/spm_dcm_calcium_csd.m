function DCM = spm_dcm_calcium_csd(P)
%==========================================================================
% Estimate parameters of a DCM using cross spectral densities (Calcium)
% FORMAT DCM = spm_dcm_calcium_csd(DCM)
%
% Calcium-Imaging DCM (spectral domain)
%--------------------------------------------------------------------------
% This routine estimates effective connectivity among regions based on
% calcium fluorescence signals, using cross-spectral densities and
% Variational Laplace under dynamic calcium equations.
%
% Example:
% subject_id = 1;
% DCM = struct();
% n_nodes = size(Y,2);
% A=0.001*ones(n_nodes,n_nodes);
% DCM.a = A;
% DCM.Y.y  = Y;
% DCM.Y.dt = TR;
% DCM.name = sprintf('DCM_cai_subject_%d', subject_id);
% DCM_est = spm_dcm_calcium_csd(DCM);

% See also: spm_dcm_fmri_csd, spm_fx_calcium, spm_gx_calcium
%__________________________________________________________________________
%
% Park HJ, 2025
%==========================================================================

% Load DCM structure
%--------------------------------------------------------------------------
if isstruct(P)
    DCM = P;
    try, DCM.name; catch, DCM.name = sprintf('DCM_Ca_%s',date); end
    P   = DCM.name;
else
    load(P)
end

% Check options and specification
%--------------------------------------------------------------------------
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.analysis;   catch, DCM.options.analysis   = 'CSD'; end
try, DCM.options.order;      catch, DCM.options.order      = 4;     end
try, DCM.options.nograph;    catch, DCM.options.nograph    = spm_get_defaults('cmdline');  end
try, DCM.options.maxit;      catch, DCM.options.maxit      = 128;   end
try, DCM.options.maxnodes;   catch, DCM.options.maxnodes   = 12;     end
try, DCM.options.modality;   catch, DCM.options.modality   = 'Ca';  end
try, DCM.options.Fdcm;       catch, DCM.options.Fdcm       = [0.01 0.2]; end
try, DCM.options.Nmax;       catch, DCM.options.Nmax       = 32;    end

% Sizes and data fields
%--------------------------------------------------------------------------
try, DCM.U.u; catch, DCM.U.u = [];            end
try, DCM.n;   catch, DCM.n = size(DCM.a,1);   end
try, DCM.v;   catch, DCM.v = size(DCM.Y.y,1); end

% Analysis options
%--------------------------------------------------------------------------
DCM.options.induced    = 1;
DCM.options.nonlinear  = 0;
DCM.options.stochastic = 0;

% Detrend and scale data
%--------------------------------------------------------------------------
DCM.Y.y = spm_detrend(DCM.Y.y);
if DCM.options.centre
    DCM.U.u = spm_detrend(DCM.U.u);
end

vY          = spm_vec(DCM.Y.y);
scale       = 1/std(vY);
DCM.Y.y     = spm_unvec(vY*scale,DCM.Y.y);
DCM.Y.scale = scale;

% Handle multi-session and missing inputs
%--------------------------------------------------------------------------
n = DCM.n;
if iscell(DCM.Y.y)
    if size(DCM.b,3) ~= numel(DCM.Y.y)
        for i = 1:numel(DCM.Y.y)
            DCM.b(:,:,i)  = DCM.b(:,:,1);
        end
    end
else
    DCM.b = zeros(n,n,0);
end
if ~isfield(DCM,'c') || isempty(DCM.c) || isempty(DCM.U.u)
    DCM.c      = zeros(n,0);
    DCM.U.u    = zeros(DCM.v,1);
    DCM.U.name = {'null'};
end
DCM.d = zeros(n,n,0);

% Priors (and initial states)
%--------------------------------------------------------------------------
[pE,pC,x] = spm_dcm_calcium_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);

hE = 8;
hC = 1/128;

% Pre-specified priors (optional)
try, pE = DCM.options.pE; end
try, pC = DCM.options.pC; end
try, hE = DCM.options.hE; end
try, hC = DCM.options.hC; end

% Eigenvector reduction for large models
if n > DCM.options.maxnodes
    try
        y = DCM.Y.y - DCM.Y.X0*(pinv(DCM.Y.X0)*DCM.Y.y);
    catch
        y = spm_detrend(DCM.Y.y);
    end
    V = spm_svd(y');
    V = V(:,1:DCM.options.maxnodes);
    j = 1:(n*n);
    V = kron(V*V',V*V');
    % Convert sparse pC to full for indexed assignment, then back to sparse
    pC_vec = spm_vec(pC);
    pC_vec(j) = V*pC_vec(j);
    pC = spm_unvec(pC_vec, pC);
end

% place eigenmodes in model if DCM.a is a vector (of eigenvalues)
%--------------------------------------------------------------------------
if isvector(DCM.a)
    DCM.M.modes = spm_svd(cov(DCM.Y.y));
end

% Model specification
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_calcium_mtf';
DCM.M.g  = @spm_gx_calcium;   % default Hill obs
if isfield(DCM.options,'linear_obs') && DCM.options.linear_obs
    DCM.M.g = @(x,u,P,M) x(:,2);  % linear: Ca concentration (second state)
end
if isfield(DCM.options,'custom_g')
    DCM.M.g = DCM.options.custom_g;
end

% Allow user to specify custom state function (e.g., spm_fx_calcium_fast)
try
    if ~isfield(DCM.M, 'f') || isempty(DCM.M.f)
        DCM.M.f  = @spm_fx_calcium;  % Default
    end
    % If user provided string, convert to handle
    if ischar(DCM.M.f) || isstring(DCM.M.f)
        DCM.M.f = str2func(DCM.M.f);
    end
catch
    DCM.M.f  = @spm_fx_calcium;  % Fallback to default
end

DCM.M.x  = x;
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = hE;
DCM.M.hC = hC;
DCM.M.n  = length(spm_vec(x));
DCM.M.m  = size(DCM.U.u,2);
DCM.M.l  = n;
DCM.M.p  = DCM.options.order;
DCM.M.nograph = DCM.options.nograph;
DCM.M.noprint = ~spm_get_defaults('dcm.verbose');

% specify M.u - endogenous fluctuations
DCM.M.u = sparse(n,1);

% Compute cross-spectral features (MAR model)
%--------------------------------------------------------------------------
DCM.Y.p  = DCM.M.p;
DCM      = spm_dcm_calcium_csd_data(DCM);  % same routine, works generically
DCM.M.Hz = DCM.Y.Hz;
DCM.M.dt = DCM.Y.dt;
DCM.M.N = length(DCM.Y.Hz);
DCM.M.ns = 1/DCM.Y.dt;

% Spectral observation noise precision
%--------------------------------------------------------------------------
DCM.Y.Q  = spm_dcm_csd_Q(DCM.Y.csd);
DCM.Y.X0 = sparse(size(DCM.Y.Q,1),0);

% Variational Laplace Inversion
%--------------------------------------------------------------------------
Y        = DCM.Y;
Y.y      = DCM.Y.csd;
[Ep,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.U,Y);

% Posterior inference
%--------------------------------------------------------------------------
dp = spm_vec(Ep) - spm_vec(pE);
d = diag(Cp);
d(d < 1e-10) = 1e-10;
Cp(1:size(Cp,1)+1:end) = d;
sd = sqrt(max(diag(Cp), 1e-16));
Pp = spm_unvec(1 - spm_Ncdf(0, abs(dp), sd), Ep);


% Predictions and residuals (capture direct outputs for debugging)
%--------------------------------------------------------------------------
% spm_csd_calcium_mtf returns [y,w,S,Gu,Gn]
[Hc, Hz_out, S_direct, Gu_direct, Gn_direct] = spm_csd_calcium_mtf(Ep,DCM.M,DCM.U);
Ec = spm_unvec(spm_vec(Y.y) - spm_vec(Hc),Hc);

% Directed transfer functions
M      = DCM.M;
Qp     = Ep;
Qp.B   = zeros(n,n,0);
Qp.C   = speye(n,n);
[S,H1] = spm_dcm_mtf_calcium(Qp,M);

% Neuronal kernel
M.g    = @(x,u,P,M) x(:,1);
[S,K1] = spm_dcm_mtf_calcium(Qp,M);

% Spectral features
% predictions (at the level of neuronal interactions)
%--------------------------------------------------------------------------
%Qp.b        = Qp.b - 32;                             % Switch off noise
%Qp.c        = Qp.c - 32;                             % Switch off noise
Qp.C        = Ep.C;
% Capture Qp-based spectral predictions and diagnostics
[Hs,Hz,dtf,Gu_qp,Gn_qp] = spm_csd_calcium_mtf(Qp,M,DCM.U);
[ccf,pst]   = spm_csd2ccf(Hs,Hz);
[coh,fsd]   = spm_csd2coh(Hs,Hz);

% Store results (include direct transfer/spectral outputs for debugging)
%--------------------------------------------------------------------------
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.Pp = Pp;
DCM.Hc = Hc;                 % Predicted (processed) CSD used in inversion
DCM.Hz_out = Hz_out;         % Frequency vector from spm_csd_calcium_mtf
DCM.S_direct = S_direct;     % Directed transfer functions (complex)
DCM.Gu_direct = Gu_direct;   % Neuronal fluctuation spectra used
DCM.Gn_direct = Gn_direct;   % Observation noise spectra used
DCM.Rc = Ec;                 % Residuals (observed - predicted)
DCM.Hs = Hs;                 % Qp-based spectral prediction
DCM.Gu_qp = Gu_qp;           % Qp neuronal spectra
DCM.Gn_qp = Gn_qp;           % Qp noise spectra
DCM.Ce = exp(-Eh);
DCM.F  = F;
DCM.dtf = dtf;
DCM.ccf = ccf;
DCM.coh = coh;
DCM.fsd = fsd;
DCM.pst = pst;
DCM.Hz  = Hz;
DCM.H1  = H1;
DCM.K1  = K1;

% Save version and results
%--------------------------------------------------------------------------
[DCM.version.SPM.version, DCM.version.SPM.revision] = spm('Ver');
DCM.version.DCM.version  = spm_dcm_ui('Version');
DCM.version.DCM.revision = 'CaI';

save(P,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
end
