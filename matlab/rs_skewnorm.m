function r = rs_skewnorm(S,ws,us,sigmas,lambdas)
%rs posterior probability of xi given S and parameters
%   Detailed explanation goes here
    p = func_likelihood_skewnorm_mixture(S,ws,us,sigmas,lambdas);
    
    M = numel(ws);
    n = size(S, 2);
    ndims = numel(us);
    assert(ndims == numel(sigmas));
    assert(ndims == numel(lambdas));
    
    p_0 = sum(p,1);
    r = zeros(M,n);
    for j = 1:M
        r(j,:) = p(j,:) ./ p_0;
    end

end

