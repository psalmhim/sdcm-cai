function y = spm_gx_linear_ab(x, u, P, M)
% Linear observation with per-region gain and offset.
%   y_i = a_i * Ca_i + b_i
%
% P.a_obs  [N x 1] log-scaled regional gain  (actual a_i = exp(P.a_obs_i))
% P.b_obs  [N x 1] linear regional offset
%
% At prior mean (pE.a_obs=0, pE.b_obs=0): y = Ca  (identical to linear_obs)

Ca = x(:, 2);                    % Ca state, all nodes [N x 1]
a  = exp(P.a_obs(:));            % regional gain, always positive
b  = P.b_obs(:);                 % regional offset
y  = a .* Ca + b;
end
