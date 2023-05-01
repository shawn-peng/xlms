function p = rr(S,ws,us,sigmas,lambdas)
%rs posterior probability of xi given s and parameters
%   Detailed explanation goes here
    assert(sum(ws) - 1.0 <= 1e-4);
    n = size(ws,2);
    assert(size(us,2) == n);
    assert(size(sigmas,2) == n);
    assert(size(lambdas,2) == n);
    
    m = size(s,2);
    p_1 = zeros(n,m);
    for i = 1:n
        w = ws(i);
        u = us(i);
        sigma = sigmas(i);
        lambda = lambdas(i);
        p_1(i,:) = w * skew_norm_pdf(s, u, sigma, lambda);
    end
    p_0 = sum(p_1,1);
    for i = 1:n
        p(i,:) = p_1(i,:) ./ p_0;
    end
end

