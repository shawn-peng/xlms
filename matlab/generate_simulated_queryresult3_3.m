function mat = generate_simulated_queryresult3_3(M, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3)
%generate_simulated_queryresult Generate a simulated dataset
%   
n = M;
% m = N;
c = randn_skew([n, 1], u_c, sigma_c, lambda_c);
i1 = randn_skew([n, 1], u_i1, sigma_i1, lambda_i1);
i2 = randn_skew([n, 1], u_i2, sigma_i2, lambda_i2);
i3 = randn_skew([n, 1], u_i3, sigma_i3, lambda_i3);

z1 = rand(n, 1) < alpha;
z2 = rand(n, 1) < beta;

s1 = z1 .* c + (1 - z1) .* i1;
s2 = z1 .* i1 + (1 - z1) .* z2 .* c + (1 - z1) .* (1 - z2) .* i2;
s3 = (z1 + (1 - z1) .* z2) .* i2 + (1 - z1) .* (1 - z2) .* i3;

mat = [s1, s2, s3];

end

