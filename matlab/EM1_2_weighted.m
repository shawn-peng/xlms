function [w1, w2, u_1, sigma_1, lambda_1, u_2, sigma_2, lambda_2] = EM1_2_weighted(s, weights, ws, u_1, sigma_1, lambda_1, u_2, sigma_2, lambda_2)

tolerance = 1e-4;

n = size(s, 2);
sum_weights = sum(weights);
proportion = sum_weights / n;

dll = inf;
ll = func_ll_skewnorm_mixture_weighted(s, weights, ws, ...
    [u_1, u_2], ...
    [sigma_1, sigma_2], ...
    [lambda_1, lambda_2]);

fig_fit = figure(1001);
ax_fit = gca;
cla(ax_fit);
hold on;
while dll > tolerance
    prev_ll = ll;
    Rs = rs_skewnorm(s, ws, ...
        [u_1, u_2], ...
        [sigma_1, sigma_2], ...
        [lambda_1, lambda_2]);
    
    Rs1 = Rs(1,:) .* weights;
    Rs2 = Rs(2,:) .* weights;
    
    cla;
    weighted_hist(s, weights, 100)
    weighted_hist(s, Rs1, 100)
    plot_skewnorm(ax_fit, s, u_1, sigma_1, lambda_1, ws(1))
    plot_skewnorm(ax_fit, s, u_2, sigma_2, lambda_2, ws(2))
    plot_skewnorm(ax_fit, s, ...
        [u_1, u_2], ...
        [sigma_1, sigma_2], ...
        [lambda_1, lambda_2], ws)
    pause(0.005);
    
    w1 = sum(Rs1) / n;
    w2 = sum(Rs2) / n;
    ws = [w1, w2];
    
    [u_1, sigma_1, lambda_1] = fit_skewnorm_weighted_one_step(s, Rs1, u_1, sigma_1, lambda_1);
    [u_2, sigma_2, lambda_2] = fit_skewnorm_weighted_one_step(s, Rs2, u_2, sigma_2, lambda_2);
    
    ll = func_ll_skewnorm_mixture_weighted(s, weights, ws, ...
        [u_1, u_2], ...
        [sigma_1, sigma_2], ...
        [lambda_1, lambda_2]);
    dll = ll - prev_ll
end

end