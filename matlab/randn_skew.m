function V = randn_skew(sz,u,sigma,lambda)
%randn_skew generate a random matrix follow a skewed normal distribution
%   Detailed explanation goes here
delta = lambda / sqrt(1+lambda^2);
Delta = sigma * delta;
Gamma = sigma^2 * (1 - delta^2);
M = randn(sz);
U = randn(sz);
V = u + Delta * abs(M) + sqrt(Gamma) * U;
end

