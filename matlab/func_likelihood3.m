
function p = func_likelihood3(S1, S2, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2)

    
    P_c1 = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_i1 = skew_norm_pdf(S1, u_i1, sigma_i1, lambda_i1);
    P_c2 = skew_norm_pdf(S2, u_c, sigma_c, lambda_c);
    P_i21 = skew_norm_pdf(S2, u_i1, sigma_i1, lambda_i1);
    P_i22 = skew_norm_pdf(S2, u_i2, sigma_i2, lambda_i2);
    
    p1 = alpha * P_c1 + (1-alpha) * P_i1;
    p2 = beta * P_c2 + alpha * P_i21 + (1-alpha-beta) * P_i22;
    
    p = [p1; p2];
end
