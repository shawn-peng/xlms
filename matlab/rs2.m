function p = rs2(s,alpha,beta,u_i1,sigma_i1,lambda_i1,u_i2,sigma_i2,lambda_i2,u_c,sigma_c,lambda_c)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    p_c = skew_norm_pdf(s, u_c, sigma_c, lambda_c);
    p_i1 = skew_norm_pdf(s, u_i1, sigma_i1, lambda_i1);
    p_i2 = skew_norm_pdf(s, u_i2, sigma_i2, lambda_i2);
    p_1 = beta * p_c;
    p_0 = beta * p_c + (1 - alpha) * p_i1 + (1 - alpha - beta) * p_i2;
    p = p_1 ./ p_0;
end

