function fdr = fdr_a(alpha, u_c, sigma_c, a_i, b_i, gamma_i, x)
%fdr_x Summary of this function goes here
%   Detailed explanation goes here

tp = alpha * (1 - normcdf(x, u_c, sigma_c));
fp = (1 - alpha) * (1 - gamcdf(x-gamma_i, a_i, 1/b_i));

fdr = fp ./ (tp + fp);


end

