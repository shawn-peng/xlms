
function ll = func_ll2_1(S1, alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1)

    
    P_c1 = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_i1 = skew_norm_pdf(S1, u_i1, sigma_i1, lambda_i1);
    
    p1 = alpha * P_c1 + (1-alpha) * P_i1;
    
    p1(p1==0) = min(p1(p1~=0));
    ll = mean(log(p1));
end
