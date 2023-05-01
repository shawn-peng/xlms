function mode = skew_norm_mode(u, sigma, lambda)

delta = lambda / sqrt(1 + lambda^2);
u_z = sqrt(2/pi) * delta;
sigma_z = sqrt(1 - u_z^2);
gamma_1 = (4-pi)/2 * (delta*sqrt(2/pi))^3 / (1 - 2*delta^2/pi)^(3/2);
m0 = u_z - gamma_1 * sigma_z / 2 + sign(lambda)/2 * exp(-2*pi/abs(lambda));
mode = u + sigma * m0;
end