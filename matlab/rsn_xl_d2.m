function r = rsn_xl_d2(s1,s2,alpha,beta,us,sigmas,lambdas)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
%   p_cic + p_ci1 + p_ici1 + p_i1i2
    M = 4;
    assert(size(us,2) == M);
    assert(size(sigmas,2) == M);
    assert(size(lambdas,2) == M);
    
    n = size(s1,2);
    assert(size(s2,2) == n);
    
    p = zeros(M,n);
    
    [u_c, u_ic, u_i1, u_i2] = unpack4(us);
    [sigma_c, sigma_ic, sigma_i1, sigma_i2] = unpack4(sigmas);
    [lambda_c, lambda_ic, lambda_i1, lambda_i2] = unpack4(lambdas);
    
    p_c_s1  = skew_norm_pdf(s1, u_c,  sigma_c,  lambda_c);
    p_ic_s1 = skew_norm_pdf(s1, u_ic, sigma_ic, lambda_ic);
    p_i1_s1 = skew_norm_pdf(s1, u_i1, sigma_i1, lambda_i1);
    
    p_c_s2  = skew_norm_pdf(s2, u_c,  sigma_c,  lambda_c);
    p_ic_s2 = skew_norm_pdf(s2, u_ic, sigma_ic, lambda_ic);
    p_i1_s2 = skew_norm_pdf(s2, u_i1, sigma_i1, lambda_i1);
    p_i2_s2 = skew_norm_pdf(s2, u_i2, sigma_i2, lambda_i2);
    
    p_cic  = alpha     * beta     * (p_c_s1  .* p_ic_s2 + p_ic_s1 .* p_c_s2 );
    p_ci1  = alpha     * (1-beta) * (p_c_s1  .* p_i1_s2 + p_i1_s1 .* p_c_s2 );
    p_ici1 = (1-alpha) * beta     * (p_ic_s1 .* p_i1_s2 + p_i1_s1 .* p_ic_s2);
	p_i1i2 = (1-alpha) * (1-beta) * (p_i1_s1 .* p_i2_s2);

	p_s = p_cic + p_ci1 + p_ici1 + p_i1i2;
	r(1,:) = p_cic  ./ p_s;
	r(2,:) = p_ci1  ./ p_s;
	r(3,:) = p_ici1 ./ p_s;
	r(4,:) = p_i1i2 ./ p_s;

end

