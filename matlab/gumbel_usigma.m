function [u, sigma] = gumbel_usigma(a, b)

c = double(eulergamma);
u = a + c * b;
sigma = b * pi / sqrt(6);

end