function p = rs1_3(s,alpha,beta,u_0,sigma_0,lambda_0,u_1,sigma_1,lambda_1,u_2,sigma_2,lambda_2)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    p_0 = skew_norm_pdf(s, u_0, sigma_0, lambda_0);
    p_1 = skew_norm_pdf(s, u_1, sigma_1, lambda_1);
    p_2 = skew_norm_pdf(s, u_2, sigma_2, lambda_2);
    p = [alpha * p_0; beta * p_1; (1-alpha-beta)*p_2];
    p_all = sum(p, 1);
    p = p ./ p_all;
end

