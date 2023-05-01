
function ll = func_ll2_1a(S1, alpha, u_c, sigma_c, a_i1, b_i1, d_i1)

    
    P_c1 = normpdf(S1, u_c, sigma_c);
    P_i1 = gampdf(S1-d_i1, a_i1, 1/b_i1);
    
    p1 = alpha * P_c1 + (1-alpha) * P_i1;
    
    p1(p1==0) = min(p1(p1~=0));
    ll = mean(log(p1));
end
