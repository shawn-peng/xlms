
function ll = func_ll(S, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i)

    [N, M] = size(S);
%     assert(N==2);
    S1 = S(1,:);
    S2 = S(2,:);
    
    P_c = skew_norm_pdf(S, u_c, sigma_c, lambda_c);
    P_i = skew_norm_pdf(S, u_i, sigma_i, lambda_i);
    
    p1 = alpha * P_c(1,:) + (1-alpha) * P_i(1,:);
    p2 = beta * P_c(2,:) + (1-beta) * P_i(2,:);
    
    p1(p1==0) = min(p1(p1~=0));
    p2(p2==0) = min(p2(p2~=0));
    ll = sum(log(p1)) + sum(log(p2));
end
