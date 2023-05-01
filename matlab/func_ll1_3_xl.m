
function ll = func_ll1_3_xl(S1, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_ic, sigma_ic, lambda_ic)

    
    P_c = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_i = skew_norm_pdf(S1, u_i, sigma_i, lambda_i);
    P_ic = skew_norm_pdf(S1, u_ic, sigma_ic, lambda_ic);
    
    p1 = alpha * P_c + beta * P_i + (1-alpha-beta) * P_ic;
    
    p1(p1==0) = min(p1(p1~=0));
    ll = mean(log(p1));
end
