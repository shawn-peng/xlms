
function ll = func_ll_skewnorm_mixture(S, ws, us, sigmas, lambdas)
% Log-likelihood of n-dimensional independent skew-normal mixtures
    p = func_likelihood_skewnorm_mixture(S, ws, us, sigmas, lambdas);
    ll = mean(log(sum(p, 1)));
end

