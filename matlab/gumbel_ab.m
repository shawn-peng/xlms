function [a, b] = gumbel_ab(u, sigma)

c = double(eulergamma);
b = sigma * sqrt(6) / pi;
a = u - c * b;

end