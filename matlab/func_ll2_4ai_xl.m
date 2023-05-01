
function ll = func_ll2_4ai_xl(S1, S2, ws, ...
    u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2)

    
    P_c1 = normpdf(S1, u_c, sigma_c);
    P_ic1 = normpdf(S1, u_ic, sigma_ic);
    P_i11 = normpdf(S1, u_i1, sigma_i1);
    P_c2 = normpdf(S2, u_c, sigma_c);
    P_ic2 = normpdf(S2, u_ic, sigma_ic);
    P_i12 = normpdf(S2, u_i1, sigma_i1);
    P_i22 = normpdf(S2, u_i2, sigma_i2);
    
    p1 = (ws(1)) * P_c1 ...
       + (ws(2)) * P_ic1 ...
       + (ws(3)) * P_i11;
    p2 = (ws(4)) * P_c2 ...
       + (ws(5)) * P_ic2 ...
       + (ws(6)) * P_i12 ...
       + (ws(7)) * P_i22;
    
    p1(p1==0) = min(p1(p1~=0));
    p2(p2==0) = min(p2(p2~=0));
    ll = mean(log(p1)) + mean(log(p2));
end
