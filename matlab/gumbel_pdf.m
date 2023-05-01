function p = gumbel_pdf(x, u, sigma)
%skew_norm_pdf the probability density at x
%   Detailed explanation goes here
[a, b] = gumbel_ab(u, sigma);
s = -(x - a) / b;
p = 1/b * exp(s - exp(s));
end

