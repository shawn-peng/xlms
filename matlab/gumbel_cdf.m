function p = gumbel_cdf(x, u, sigma)
%gumbel_cdf the cumulative density at x
%   Detailed explanation goes here
[a, b] = gumbel_ab(u, sigma);
s = -(x - a) / b;
p = exp(-exp(s));
end

