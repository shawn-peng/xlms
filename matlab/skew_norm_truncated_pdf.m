function pt = skew_norm_truncated_pdf(x, u, sigma, lambda, t)
%skew_norm_pdf the probability density at x
%   Detailed explanation goes here
flags = x <= t;

p = (2/sigma) * normpdf((x-u) / sigma) .* normcdf(lambda * (x-u) / sigma);

pt = flags .* p ./ skew_norm_cdf(t, u, sigma, lambda);

end

