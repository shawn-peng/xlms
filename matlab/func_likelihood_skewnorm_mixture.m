
function p = func_likelihood_skewnorm_mixture(S, ws, us, sigmas, lambdas)
% likelihood of n-dimensional independent skew-normal mixtures
    ndims = size(S, 1);
    assert(ndims == size(us,1));
    assert(ndims == size(sigmas,1));
    assert(ndims == size(lambdas,1));
    
    M = size(ws, 2);
    assert(size(us, 2) == M);
    assert(size(sigmas, 2) == M);
    assert(size(lambdas, 2) == M);
    
    n = size(S, 2);
    
%     lower_bond = 1e-27;
    lower_bond = 0;
    
    pp = zeros(ndims,M,n);
    p = zeros(M,n);
    for j = 1:M
        w = ws(j);
        u = us(:,j);
        sigma = sigmas(:,j);
        lambda = lambdas(:,j);
        
        for k = 1:ndims
            pp(k,j,:) = skew_norm_pdf(S(k, :), u(k), sigma(k), lambda(k));
        end
        p(j,:) = w * prod(pp(:,j,:), 1); % pp reduced on the first dim
%         p(j, (p(j,:) < lower_bond)) = lower_bond;
    end

end

