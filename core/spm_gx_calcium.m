function [y, dgdx, dgdP] = spm_gx_calcium(x, u, P, M)
%==========================================================================
% Observation function for Calcium Imaging DCM
% Maps intracellular [Ca2+] → fluorescence signal (ΔF/F)
%
% FORMAT [y, dgdx, dgdP] = spm_gx_calcium(x, u, P, M)
%
% INPUT:
%   x - state vector [neuronal; calcium] (n×2 or 2n×1)
%   u - input (unused)
%   P - parameters (alpha, Kd, n, beta_y)
%   M - model structure
%
% OUTPUT:
%   y    - fluorescence observations (n×1)
%   dgdx - Jacobian w.r.t. states
%   dgdP - Jacobian w.r.t. parameters
%
% OBSERVATION MODEL:
%   F(Ca) = Ca^n / (Kd^n + Ca^n)    [Hill equation]
%   y = α * F(Ca) + β                [fluorescence]
%
% where:
%   α     - fluorescence sensitivity (baseline=1.0)
%   Kd    - half-saturation [Ca2+] (baseline=0.5)
%   n     - Hill coefficient (baseline=2.0)
%   β     - baseline fluorescence offset
%
%--------------------------------------------------------------------------
%   Park & ChatGPT (2025)
%==========================================================================

if nargin < 4, M = struct(); end
if size(x,2) < 2, x(:,2) = zeros(size(x,1),1); end
n = size(x,1);

%--------------------------------------------------------------------------
% Extract calcium state
Ca = max(min(x(:,2), 10), 1e-6);  % prevent extreme values

%--------------------------------------------------------------------------
% Parameter exponentiation (positive parameters)
% alpha_base calibrated to match M1 linear obs gain (dg/dCa=1 at Ca_eq=0.12, Kd=0.2)
% alpha_base = 1/Hill_slope_at_eq = 1/(n*Kd/(Kd+Ca)^2) = 1/1.953 ≈ 0.512
H = [0.512 0.2 1.0];    % base scaling for [alpha, Kd, nH]
alpha = H(1) * exp(P.alpha(:));  
Kd    = H(2) * exp(P.Kd(:));    
nH    = H(3) * exp(P.n(:));     
beta  = P.beta_y; 

%--------------------------------------------------------------------------
% Observation model: Hill-type fluorescence response
F = (Ca.^nH) ./ (Kd.^nH + Ca.^nH);
y = alpha .* F + beta;

%==========================================================================
% (I) Jacobian wrt states (dgdx)
%==========================================================================
if nargout > 1
    % dF/dCa (derivative of Hill equation)
    dFdCa = (nH .* Kd.^nH .* Ca.^(nH - 1)) ./ (Kd.^nH + Ca.^nH).^2;
    
    % ∂y/∂x1 = 0 (neuronal state doesn't affect fluorescence directly)
    % ∂y/∂x2 = α * dF/dCa (calcium state affects fluorescence via Hill eq)
    dgdx = [sparse(n,n), diag(alpha .* dFdCa)];
else
    dgdx = [];
end

%==========================================================================
% (II) Jacobian wrt parameters (dgdP)
%==========================================================================
if nargout > 2
    dgdP = struct();
    
    % ∂y/∂alpha = F (fluorescence response)
    dgdP.alpha = diag(F);
    
    % ∂y/∂beta = 1 (baseline offset)
    dgdP.beta  = speye(n);
    
    % ∂y/∂Kd (half-saturation parameter)
    dFdKd = -nH .* Kd.^(nH - 1) .* (Ca.^nH) ./ (Kd.^nH + Ca.^nH).^2;
    dgdP.Kd = diag(alpha .* dFdKd);
    
    % ∂y/∂n (Hill coefficient)
    dFdN = F .* (1 - F) .* log(Ca ./ Kd);
    dgdP.n = diag(alpha .* dFdN);
else
    dgdP = [];
end
end
