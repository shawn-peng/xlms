
function ll = func_ll2_2(S1, S2, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2)

    
    P_c1 = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_i1 = skew_norm_pdf(S1, u_i1, sigma_i1, lambda_i1);
    P_c2 = skew_norm_pdf(S2, u_c, sigma_c, lambda_c);
    P_i2 = skew_norm_pdf(S2, u_i2, sigma_i2, lambda_i2);
    
    p1 = alpha * P_c1 + (1-alpha) * P_i1;
    p2 = beta * P_c2 + (1-beta) * P_i2;
    
    p1(p1==0) = min(p1(p1~=0));
    p2(p2==0) = min(p2(p2~=0));
    ll = sum(log(p1)) + sum(log(p2));
end
