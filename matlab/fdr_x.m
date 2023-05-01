function fdr = fdr_x(alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, x)
%fdr_x Summary of this function goes here
%   Detailed explanation goes here

tp = alpha * (1 - skew_norm_cdf(x, u_c, sigma_c, lambda_c));
fp = (1 - alpha) * (1 - skew_norm_cdf(x, u_i, sigma_i, lambda_i));

fdr = fp ./ (tp + fp);


end

