function [ws, us, sigmas, lambdas] = EM1_n_weighted(s, weights, ws, us, sigmas, lambdas)

tolerance = 1e-4;

n = size(s, 2);
M = numel(ws);
assert(M == numel(us));
assert(M == numel(sigmas));
assert(M == numel(lambdas));

sum_weights = sum(weights);
% proportion = sum_weights / n;

dll = inf;
ll = func_ll_skewnorm_mixture_weighted(s, weights, ws, ...
    us, sigmas, lambdas);

fig_fit = figure(1002);
ax_fit = gca;
cla(ax_fit);
hold on;
while dll > tolerance
    prev_ll = ll;
    Rs = rs_skewnorm(s, ws, us, sigmas, lambdas);
    
    cla;
    weighted_hist(s, weights, 100)
    weighted_hist(s, Rs(1,:), 100)
    for j = 1:M
        plot_skewnorm(ax_fit, s, us(j), sigmas(j), lambdas(j), ws(j))
    end
    plot_skewnorm(ax_fit, s, us, sigmas, lambdas, ws)
    pause(0.005);
    
    ws = sum(Rs, 2)' / n;
    
    for j = 1:M
        [uj, sigmaj, lambdaj] = fit_skewnorm_weighted_one_step(s, Rs(j,:), us(j), sigmas(j), lambdas(j));
        us(j) = uj;
        sigmas(j) = sigmaj;
        lambdas(j) = lambdaj;
    end
    
    ll = func_ll_skewnorm_mixture_weighted(...
        s, weights, ws, us, sigmas, lambdas);
    dll = ll - prev_ll
end

end