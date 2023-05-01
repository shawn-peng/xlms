
function ll = func_ll3_5(S1, S2, S3, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3)

    
    P_c1 = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_i1 = skew_norm_pdf(S1, u_i1, sigma_i1, lambda_i1);
    P_c2 = skew_norm_pdf(S2, u_c, sigma_c, lambda_c);
    P_i21 = skew_norm_pdf(S2, u_i1, sigma_i1, lambda_i1);
    P_i22 = skew_norm_pdf(S2, u_i2, sigma_i2, lambda_i2);
    P_i32 = skew_norm_pdf(S3, u_i2, sigma_i2, lambda_i2);
    P_i33 = skew_norm_pdf(S3, u_i3, sigma_i3, lambda_i3);
    
    p1 = alpha * P_c1 + (1-alpha) * P_i1;
    p2 = beta * P_c2 + alpha * P_i21 + (1-alpha-beta) * P_i22;
    p3 = (alpha+beta) * P_i32 + (1-alpha-beta) * P_i33;
    
    p1(p1==0) = min(p1(p1~=0));
    p2(p2==0) = min(p2(p2~=0));
    p3(p3==0) = min(p3(p3~=0));
    ll = sum(log(p1)) + sum(log(p2)) + sum(log(p3));
end
