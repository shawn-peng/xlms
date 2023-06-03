function fdr = fdr_xl(x, ...
    w_c, w_ic, w_i1, ...
    u_c, sigma_c, lambda_c, ...
    u_i, sigma_i, lambda_i, ...
    u_ic, sigma_ic, lambda_ic)
%fdr_x Summary of this function goes here
%   Detailed explanation goes here

tp = w_c * (1 - skew_norm_cdf(x, u_c, sigma_c, lambda_c));
fpic = w_ic * (1 - skew_norm_cdf(x, u_ic, sigma_ic, lambda_ic));
fpi = w_i1 * (1 - skew_norm_cdf(x, u_i, sigma_i, lambda_i));
fp = fpi + fpic;

fdr = fp ./ (tp + fp);


end

