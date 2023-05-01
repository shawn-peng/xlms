function r = rsn_xl_d2i(s1,s2,ws,us,sigmas,lambdas)
%rs posterior probability of xi given S and parameters
%   Detailed explanation goes here
    M = size(ws, 2);
    assert(size(us, 2) == M);
    assert(size(sigmas, 2) == M);
    assert(size(lambdas, 2) == M);
    
    n = size(s1, 2);
    assert(n == size(s2, 2));
    
%     lower_bond = 1e-27;
    lower_bond = 0;
    
    p1 = zeros(M,n);
    p2 = zeros(M,n);
    p = zeros(M,n);
    for j = 1:M
        w = ws(j);
        u = us(:,j);
        sigma = sigmas(:,j);
        lambda = lambdas(:,j);
        
        p1(j,:) = skew_norm_pdf(s1, u(1), sigma(1), lambda(1));
        p2(j,:) = skew_norm_pdf(s2, u(2), sigma(2), lambda(2));
        p(j,:) = w * p1(j,:) .* p2(j,:);
%         p(j, (p(j,:) < lower_bond)) = lower_bond;
    end
 
    p_0 = sum(p,1);
    r = zeros(M,n);
    for j = 1:M
        r(j,:) = p(j,:) ./ p_0;
    end
end

