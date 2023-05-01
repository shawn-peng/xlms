function [alpha, beta] = fit_gamma_weighted(s, weights)

sum_w = sum(weights);

u = sum(weights .* s) / sum_w;
variance = sqrt(sum(weights .* (s - u).^2) / sum_w);

alpha = u.^2 ./ variance;
beta = u ./ variance;

end