
function ll = func_ll2_1g(S1, alpha, u_c, sigma_c, u_i1, sigma_i1)

    
    P_c1 = normpdf(S1, u_c, sigma_c);
    P_i1 = gumbel_pdf(S1, u_i1, sigma_i1);
    
    p1 = alpha * P_c1 + (1-alpha) * P_i1;
    
    p1(p1==0) = min(p1(p1~=0));
    ll = sum(log(p1));
end
