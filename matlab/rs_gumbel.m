function p = rs_gumbel(s, alpha, u_c, sigma_c, u_i, sigma_i)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    n = 2;
    
    m = size(s,2);
    p = zeros(n,m);
    p_1 = zeros(n,m);

    p_1(1,:) = alpha * normpdf(s, u_c, sigma_c);
    p_1(2,:) = (1 - alpha) * gumbel_pdf(s, u_i, sigma_i);
    
    p_0 = sum(p_1,1);
%     for i = 1:n
%         p(i,:) = p_1(i,:) ./ p_0;
%     end
    p = p_1(1,:) ./ p_0;
end