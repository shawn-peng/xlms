
Rs = rs_skewnorm(S, ws, ...
    [u_c1,  u_c1,  u_ic1, u_i11, u_ic1, u_i11, u_i11;
     u_ic2, u_i12, u_c2,  u_c2,  u_i12, u_ic2, u_i22], ...
    [sigma_c1,  sigma_c1,  sigma_ic1, sigma_i11, sigma_ic1, sigma_i11, sigma_i11;
     sigma_ic2, sigma_i12, sigma_c2,  sigma_c2,  sigma_i12, sigma_ic2, sigma_i22], ...
    [lambda_c1,  lambda_c1,  lambda_ic1, lambda_i11, lambda_ic1, lambda_i11, lambda_i11;
     lambda_ic2, lambda_i12, lambda_c2,  lambda_c2,  lambda_i12, lambda_ic2, lambda_i22]);



fig_debug_s1 = figure;
ax_debug_s1 = axes(fig_debug_s1);
hold on

Rsj1_initial = Rs1c;

% histogram(s1, 100, 'Normalization', 'pdf')
weighted_hist(s1, ones(M1,1), 100)
weighted_hist(s1, Rsj1, 100)
weighted_hist(s1, Rs1c, 100)
% weighted_hist(s1, Rsj2, 100)
% weighted_hist(s1, Rsj3, 100)
% weighted_hist(s1, Rsj4, 100)
% weighted_hist(s1, Rsj5, 100)
% weighted_hist(s1, Rsj6, 100)
weighted_hist(s1, Rsj7, 100)

% legend({'s1'; 'Rsj1'; 'Rsj2'; 'Rsj3'; 'Rsj4'; 'Rsj5'; 'Rsj6'; 'Rsj7';});

plot_skewnorm(ax_debug_s1, s1, u_c1, sigma_c1, lambda_c1, w1)

plot_skewnorm(ax_debug_s1, s1, u_ic1, sigma_ic1, lambda_ic1, w5)

[u_i11, sigma_i11, lambda_i11] = fit_skewnorm_weighted(s1, Rsj7, u_i11, sigma_i11, lambda_i11-1);
plot_skewnorm(ax_debug_s1, s1, u_i11, sigma_i11, lambda_i11, w7)

legend

fig_debug_s2 = figure;
ax_debug_s2 = axes(fig_debug_s2);
hold on

weighted_hist(s2, ones(M1,1), 100)
weighted_hist(s2, dRs1, 100)
weighted_hist(s2, rRs1, 100)
weighted_hist(s2, Rsj1, 100)
% weighted_hist(s2, Rsj1+Rsj5, 100)
weighted_hist(s2, Rsj2, 100)
weighted_hist(s2, Rsj3, 100)
weighted_hist(s2, Rsj4, 100)
weighted_hist(s2, Rsj5, 100)
weighted_hist(s2, Rsj6, 100)
weighted_hist(s2, Rsj7, 100)

legend({'s2'; 'Rsj1'; 'Rsj2'; 'Rsj3'; 'Rsj4'; 'Rsj5'; 'Rsj6'; 'Rsj7';});

plot_skewnorm(ax_debug_s2, s2, u_c2, sigma_c2, lambda_c2, w3)

% [u_ic2, sigma_ic2, lambda_ic2] = fit_skewnorm_weighted(s2, Rsj1, u_ic2, sigma_ic2, lambda_ic2)
plot_skewnorm(ax_debug_s2, s2, u_ic2, sigma_ic2, lambda_ic2, w1)

% [u_i12, sigma_i12, lambda_i12] = fit_skewnorm_weighted(s2, Rsj5, u_i12, sigma_i12, lambda_i12)
plot_skewnorm(ax_debug_s2, s2, u_i12, sigma_i12, lambda_i12, w5)
% plot_skewnorm(ax_debug_s2, s2, u_i12, sigma_i12+5, lambda_i12-1, w5)

[u_i22, sigma_i22, lambda_i22] = fit_skewnorm_weighted(s2, Rsj7, u_i22, sigma_i22, lambda_i22);
plot_skewnorm(ax_debug_s2, s2, u_i22, sigma_i22, lambda_i22-2, w7)

% [alpha, beta] = fit_gamma_weighted(s2, Rsj5)
% plot_gamma(s2, alpha, beta, w5)

figure
hold on;

scatter(s1, s2, 10, Rsj1, '.')
scatter(s1, s2, 10, Rsj2, '.')
scatter(s1, s2, 10, Rsj3, '.')
scatter(s1, s2, 10, Rsj4, '.')
scatter(s1, s2, 10, Rsj5, '.')
scatter(s1, s2, 10, Rsj6, '.')
scatter(s1, s2, 10, Rsj7, '.')
colorbar


ws(1) = (w1+w3)/2;
ws(3) = ws(1)
Rsdebug = rs_skewnorm(S, ws, ...
    [u_c1,  u_c1,  u_ic1, u_i11, u_ic1, u_i11, u_i11;
     u_ic2, u_i12, u_c1,  u_c1,  u_i12, u_ic2, u_i22], ...
    [sigma_c1,  sigma_c1,  sigma_ic1, sigma_i11, sigma_ic1, sigma_i11, sigma_i11;
     sigma_ic2, sigma_i12, sigma_c1,  sigma_c1,  sigma_i12, sigma_ic2, sigma_i22], ...
    [lambda_c1,  lambda_c1,  lambda_ic1, lambda_i11, lambda_ic1, lambda_i11, lambda_i11;
     lambda_ic2, lambda_i12, lambda_c1,  lambda_c1,  lambda_i12, lambda_ic2, lambda_i22]);

Rsj1 = Rsdebug(1,:);
Rsj2 = Rsdebug(2,:);
Rsj3 = Rsdebug(3,:);
Rsj4 = Rsdebug(4,:);
Rsj5 = Rsdebug(5,:);
Rsj6 = Rsdebug(6,:);
Rsj7 = Rsdebug(7,:);


bin80_flags = S(1,:) > 78 & S(1,:) < 82;
bin80 = S(:,bin80_flags);

figure
hold on;
weighted_hist(S(2,:), bin80_flags, 100)
weighted_hist(S(2,:), Rsj4.*bin80_flags, 100)
weighted_hist(S(2,:), Rsj5.*bin80_flags, 100)
weighted_hist(S(2,:), Rsj6.*bin80_flags, 100)
weighted_hist(S(2,:), Rsj7.*bin80_flags, 100)



