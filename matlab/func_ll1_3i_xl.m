
function ll = func_ll1_3i_xl(S1, ws, ...
    u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1)

    
    P_c1 = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_ic1 = skew_norm_pdf(S1, u_ic, sigma_ic, lambda_ic);
    P_i11 = skew_norm_pdf(S1, u_i1, sigma_i1, lambda_i1);
    
    p1 = (ws(1)) * P_c1 ...
       + (ws(2)) * P_ic1 ...
       + (ws(3)) * P_i11;
    
    p1(p1==0) = min(p1(p1~=0));
    ll = mean(log(p1));
end
