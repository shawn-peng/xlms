function p = rs_1(s,alpha,u_i,sigma_i,lambda_i,u_c,sigma_c,lambda_c)
%rs_1 posterior probability of xi given s and parameters
%   Detailed explanation goes here
    N = size(s, 1);
    p_i = skew_norm_pdf(s, u_i, sigma_i, lambda_i);
    p_c = skew_norm_pdf(s, u_c, sigma_c, lambda_c);
    p_f = p_c ./ p_i;
    p_1 = alpha * p_f;
    p_0 = alpha * sum(p_f, 1) + (1 - alpha) * N;
    p = p_1 ./ (ones(N,1) * p_0);
end

