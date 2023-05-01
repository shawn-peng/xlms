function p = rs_gamma(s, alpha, u_c, sigma_c, alpha_i, beta_i, gamma_i)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    n = 2;
    
    m = size(s,2);
    p = zeros(n,m);
    p_1 = zeros(n,m);

    p_1(1,:) = alpha * normpdf(s, u_c, sigma_c);
    p_1(2,:) = (1 - alpha) * gampdf(s-gamma_i, alpha_i, 1/beta_i);
    
    p_0 = sum(p_1,1);
%     for i = 1:n
%         p(i,:) = p_1(i,:) ./ p_0;
%     end
    p = p_1(1,:) ./ p_0;
end