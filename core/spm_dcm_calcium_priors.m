function [pE,pC,x,C] = spm_dcm_calcium_priors(A,B,C,D,options)
%==========================================================================
% Priors for calcium spectral DCM (population-level)
% Returns:
%   pE : prior expectations
%   pC : prior covariances
%   x  : initial states [neuronal, calcium]
%   C  : diagonalised prior covariance
%==========================================================================

% number of regions
n = length(A);

% options
try, pA = exp(options.precision); catch, pA = 8; end
try, options.induced; catch, options.induced = 1; end
if nargin < 4 || isempty(D), D = zeros(n,n,0); end

%--------------------------------------------------------------------------
% Initial states (consistent with v2 operating point: sigma(x1=0)=0.5)
%--------------------------------------------------------------------------
H_init  = [0.8 0.02 0.26 0.2];  % [tau, beta, kappa, R] base values
x2_ss   = H_init(1) * (H_init(3)*0.5 + H_init(2));  % SS calcium at sigma=0.5: ~0.12
x = sparse(n,2);
x(:,1) = 0.0;    % neuronal state (operating point)
x(:,2) = x2_ss;  % calcium at steady state when x1=0, sigma=0.5

%--------------------------------------------------------------------------
% Connectivity priors
%--------------------------------------------------------------------------
A = logical(A);

pE.A = A/128;
pE.A = pE.A - diag(diag(A));   % self-inhibition (log-domain)
pE.A(logical(eye(n))) = 0;     % explicit fix: self-connection prior mean always 0,
                                % independent of whatever structural mask is passed in
                                % for off-diagonal connections (previously the
                                % all-ones mask used in production also shifted the
                                % diagonal to -127/128 via the subtraction above)

pE.B = B*0;
pE.C = C*0;
pE.D = D*0;

for i = 1:size(A,3)
    pC.A(:,:,i) = double(A(:,:,i))/pA;
    pC.A(:,:,i) = pC.A(:,:,i) - diag(diag(pC.A(:,:,i))) + 2*eye(n);
    % explicit fix: self-connection prior variance always 2 (SPM default),
    % independent of the off-diagonal structural-mask-derived precision
end
pC.B = double(B)/pA;
pC.C = double(C)/pA;
pC.D = double(D)/pA;

%--------------------------------------------------------------------------
% Calcium dynamics priors
%--------------------------------------------------------------------------
pE.kappa_Ca = 0;
pC.kappa_Ca = 1/64;

pE.beta_Ca  = zeros(n,1);
pC.beta_Ca  = zeros(n,1) + 1/64;

pE.tau_Ca   = 0;
pC.tau_Ca   = 1/64;

pE.R        = 0;
pC.R        = 1/32;

pE.V0       = 0;    % v2 mapping: V0 is resting voltage (σ=0.5 at x1=0 by construction)
pC.V0       = 1/16;

%--------------------------------------------------------------------------
% Observation model (fluorescence)
%--------------------------------------------------------------------------
pE.alpha  = 0;
pC.alpha  = 1/16;

pE.beta_y = zeros(n,1);
pC.beta_y = ones(n,1) * (1/64);

pE.Kd     = 0;
pC.Kd     = 1/128;

pE.n      = 0;
pC.n      = 1/128;

%--------------------------------------------------------------------------
% Spectral priors
%--------------------------------------------------------------------------
if options.induced
    pE.a =  sparse(2,1);  
    pE.b =  sparse(2,1);  
    pE.c = sparse(1,n);

    pC.a = sparse(2,1) + 1/64;
    pC.b = sparse(2,1) + 1/64;
    pC.c = sparse(1,n) + 1/64;
end

%--------------------------------------------------------------------------
% Flatten covariance
%--------------------------------------------------------------------------
C = diag(spm_vec(pC));

end
