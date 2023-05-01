function [r1_c, r1_i1, r2_c, r2_i1, r2_i2, r3_i2, r3_i3] = rsdn(S,ws,us,sigmas,lambdas)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    assert(sum(ws) - 1.0 <= 1e-4);
%     n = size(ws,2);
%     assert(size(us,2) == n);
%     assert(size(sigmas,2) == n);
%     assert(size(lambdas,2) == n);
    
    w_1 = ws(1);
    w_2 = ws(2);
    w_3 = ws(3);
    
    u_c = us(1);
    u_i1 = us(2);
    u_i2 = us(3);
    u_i3 = us(4);
    
    sigma_c = sigmas(1);
    sigma_i1 = sigmas(2);
    sigma_i2 = sigmas(3);
    sigma_i3 = sigmas(4);
    
    lambda_c = lambdas(1);
    lambda_i1 = lambdas(2);
    lambda_i2 = lambdas(3);
    lambda_i3 = lambdas(4);
    

    s1 = S(1,:);
    s1 = s1(s1~=0);
    s2 = S(2,:);
    flags2 = s2~=0;
    s3 = S(3,:);
    flags3 = s3~=0;
    
    p1_c = w_1 * skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
    p1_i1 = (w_2 + w_3) * skew_norm_pdf(s1, u_i1, sigma_i1, lambda_i1);
    
    r1_c = p1_c ./ (p1_c + p1_i1);
    r1_i1 = p1_i1 ./ (p1_c + p1_i1);
    
    p2_c = w_1 * skew_norm_pdf(s2, u_c, sigma_c, lambda_c);
    p2_i1 = w_2 * skew_norm_pdf(s2, u_i1, sigma_i1, lambda_i1);
    p2_i2 = w_3 * skew_norm_pdf(s2, u_i2, sigma_i2, lambda_i2);
    
    r2_c = r1_i1 .* p2_c ./ (p2_c + p2_i1 + p2_i2);
    r2_i1 = r1_c .* p2_i1 ./ (p2_c + p2_i1 + p2_i2);
    r2_i2 = r1_i1 .* p2_i2 ./ (p2_c + p2_i1 + p2_i2);
    
    p3_i2 = (w_1 + w_2) * skew_norm_pdf(s3, u_i2, sigma_i2, lambda_i2);
    p3_i3 = w_3 * skew_norm_pdf(s3, u_i3, sigma_i3, lambda_i3);
    
    r3_i2 = (r2_c + r2_i1) .* p3_i2 ./ (p3_i2 + p3_i3);
    r3_i3 = r2_i2 .* p3_i3 ./ (p3_i2 + p3_i3);
    
    r2_c = r2_c(flags2);
    r2_i1 = r2_i1(flags2);
    r2_i2 = r2_i2(flags2);
    
    r3_i2 = r3_i2(flags3);
    r3_i3 = r3_i3(flags3);
end

