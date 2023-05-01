function p = skew_norm_cdf(x, u, sigma, lambda)
%skew_norm_pdf the probability density at x
%   Detailed explanation goes here
n = size(x, 1);
m = size(x, 2);
for i = 1:n
    for j = 1:m
        p(i,j) = normcdf((x(i,j)-u) / sigma) - 2 * tfn((x(i,j)-u) / sigma, lambda);
    end
end
end

