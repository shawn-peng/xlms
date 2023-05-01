function p = rs(s,alpha,u_0,sigma_0,lambda_0,u_1,sigma_1,lambda_1)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    p_1 = skew_norm_pdf(s, u_1, sigma_1, lambda_1);
    p_0 = skew_norm_pdf(s, u_0, sigma_0, lambda_0);
    p_n = alpha * p_1;
    p_all = alpha * p_1 + (1 - alpha) * p_0;
    p = p_n ./ p_all;
end

