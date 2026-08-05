function y = spm_gx_calcium_linear_ab(x, u, P, M)
% Linear-ab observation: y_i = exp(alpha) * Ca_i + beta_y_i
% Per-region offset (beta_y), global gain (alpha).
if size(x,2) < 2, x(:,2) = zeros(size(x,1),1); end
Ca    = x(:,2);
alpha = exp(P.alpha(1));
y     = alpha .* Ca + P.beta_y;
end
