function [y, dgdx, dgdP] = spm_gx_calcium_linear(x, u, P, M)
%==========================================================================
% Linear observation function for Calcium Spectral DCM (M2 in BMS v2).
%
% Replaces the Hill saturation curve with a simple linear proportionality
% to the intracellular calcium state. Tests whether the S-shaped Hill
% nonlinearity contributes beyond linear Ca-to-fluorescence scaling.
%
%   y = alpha * Ca + beta
%
% Jacobian dgdx(:,Ca) = alpha (constant slope, no saturation)
% vs Hill:  dgdx(:,Ca) = alpha * dF/dCa (slope depends on operating point)
%
% Uses same parameter struct as spm_gx_calcium.m for compatibility.
% Parameters Kd and n are ignored (Hill shape not present).
%
%==========================================================================

if nargin < 4, M = struct(); end
if size(x,2) < 2, x(:,2) = zeros(size(x,1),1); end
n = size(x, 1);

Ca    = max(min(x(:,2), 10), 1e-6);
H_alpha = 8.0;                        % same base as spm_gx_calcium.m
alpha = H_alpha * exp(P.alpha(:)) .* ones(n,1);  % broadcast scalar gain to
                                                  % n regions (P.alpha is a
                                                  % single global parameter,
                                                  % same as in spm_gx_calcium.m;
                                                  % without this broadcast,
                                                  % diag(alpha) below is 1x1
                                                  % instead of nxn, breaking
                                                  % the Jacobian dimensions)
beta  = P.beta_y;

y = alpha .* Ca + beta;

if nargout > 1
    dgdx = [sparse(n,n), diag(alpha)];   % slope = alpha, no Hill factor
else
    dgdx = [];
end

if nargout > 2
    dgdP = struct();
    dgdP.alpha  = diag(Ca);
    dgdP.beta_y = speye(n);
    dgdP.Kd     = sparse(n,n);   % unused; kept for struct compatibility
    dgdP.n      = sparse(n,n);
else
    dgdP = [];
end
end
