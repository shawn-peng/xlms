function plot_skewnorm_cdf(ax, x, u, sigma, lambda, w)
	x = min(x):0.01:max(x);
    M = numel(w);
    assert(M == numel(u));
    assert(M == numel(sigma));
    assert(M == numel(lambda));
    for i = 1:M
        y(:,i) = skew_norm_cdf(x, u(i), sigma(i), lambda(i));
    end
	plot(ax, x, y*reshape(w, [M,1]))
end
