function mat = generate_simulated_queryresult2_1a(M, alpha, u_c, sigma_c, a_i, b_i)
%generate_simulated_queryresult Generate a simulated dataset
%   
n = M;
% m = N;
c = normrnd(u_c, sigma_c, [n, 1]);
i1 = gamrnd(a_i, 1/b_i, [n, 1]);

z1 = rand(n, 1) < alpha;

s1 = z1 .* c + (1 - z1) .* i1;

mat = s1;

end

