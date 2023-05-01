function p = skew_norm_pdf(x, u, sigma, lambda)
%skew_norm_pdf the probability density at x
%   Detailed explanation goes here
p = (2/sigma) .* normpdf((x-u) ./ sigma) .* normcdf(lambda .* (x-u) ./ sigma);
end

