function [u, sigma] = fit_normal_weighted(s, weights)

sum_w = sum(weights);

u = sum(weights .* s) / sum_w;
sigma = sqrt(sum(weights .* (s - u).^2) / sum_w);

end