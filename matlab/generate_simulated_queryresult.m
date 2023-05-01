function M = generate_simulated_queryresult(N, M, alpha, u_i, sigma_i, lambda_i, u_c, sigma_c, lambda_c)
%generate_simulated_queryresult Generate a simulated dataset
%   
n = M;
m = N;
M = randn_skew([n, m], u_i, sigma_i, lambda_i);
c = randn_skew([n, 1], u_c, sigma_c, lambda_c);

cf = rand(n, 1) < alpha;
sum(cf)

ci = randi([1, m], n, 1);
ci = sub2ind([n, m], [1:n]', ci);
% ci_flag = zeros(n, m);
% ci_flag(ci) = 1;
% ci_flag = ci_flag~=0;
M(ci) = cf .* c + (1-cf) .* M(ci);
M = M';

end

