function [dlls, ws, us, sigmas, lambdas] = EM2_10j_xl(S, lambda_signs)

tolerance = 1e-8;

S = S(:, S(2,:)~=0);
q1 = quantile(S(1,:), 0.01);
S = S(:, (S(1,:) > q1));
N = size(S, 1);

prev_ll = -inf;
dll = inf;

s1 = S(1,:);
s2 = S(2,:);

M1 = size(s1, 2);
M2 = size(s2, 2);

[alpha, beta, ...
    u_c1, sigma_c1, lambda_c1, ...
    u_ic1, sigma_ic1, lambda_ic1, ...
    u_i11, sigma_i11, lambda_i11] = EM1_3i_xl(S, lambda_signs(1), lambda_signs(2), lambda_signs(3));

w1c = alpha;
w1ic = beta;
w1i1 = 1 - alpha - beta;

Rs1 = rsn(s1, ...
        [w1c, w1ic, w1i1], ...
        [u_c1, u_ic1, u_i11], ...
        [sigma_c1, sigma_ic1, sigma_i11], ...
        [lambda_c1, lambda_ic1, lambda_i11]);

Rs1c = Rs1(1, :);
Rs1ic = Rs1(2, :);
Rs1i1 = Rs1(3, :);

w12 = sum(Rs1c) / M1;
w34 = sum(Rs1ic) / M1;
w567 = sum(Rs1i1) / M1;
% s2
u_c2 = u_c1;
sigma_c2 = sigma_c1;

[u_ic2, sigma_ic2] = fit_normal_weighted(s2, Rs1c);
[u_i12, sigma_i12] = fit_normal_weighted(s2, Rs1ic);
[u_i22, sigma_i22] = fit_normal_weighted(s2, Rs1i1);

% lambda_c2  = lambda_signs(1);
% lambda_ic2 = lambda_signs(1);
% lambda_i12 = lambda_signs(3);
% lambda_i22 = lambda_signs(4);

lambda_ic2 = 1;
lambda_i12 = -1;
lambda_i22 = -1;
lambda_ic22 = 1;
lambda_i12ic = -1;
lambda_c2i1 = 1;
lambda_ic2i1 = 1;
lambda_i22 = -1;

[u_ic2c, sigma_ic2c, lambda_ic2c] = fit_skewnorm_weighted(s2, Rs1c,  u_ic2, sigma_ic2, lambda_ic2);
[u_i12ic, sigma_i12ic, lambda_i12ic] = fit_skewnorm_weighted(s2, Rs1ic, u_i12, sigma_i12, lambda_i12);
[u_i22, sigma_i22, lambda_i22] = fit_skewnorm_weighted(s2, Rs1i1, u_i22, sigma_i22, lambda_i22);

% [w1, w2, ...
%     u_c1, sigma_c1, lambda_c1, ...
%     u_ic1, sigma_ic1, lambda_ic1] ...
%     = EM1_2_weighted(s1, ones(1,M1), [w12/2, w12/2], ...
%         u_c1, sigma_c1, lambda_c1, ...
%         u_ic1, sigma_ic1, lambda_ic1);
% 
% [ws, us, sigmas, lambdas] = EM1_n_weighted(s1, ones(1,M1), ...
%     [w1c, w1ic, w1i1], ...
%     [u_c1, u_ic1, u_i11], ...
%     [sigma_c1, sigma_ic1, sigma_i11], ...
%     [lambda_c1, lambda_ic1, lambda_i11]);
    
[w1, w2, ...
    u_ic2c, sigma_ic2c, lambda_ic2c, ...
    u_i12c, sigma_i12c, lambda_i12c] ...
    = EM1_2_weighted(s2, Rs1c, [w12/2, w12/2], ...
        u_ic2, sigma_ic2, lambda_ic2, ...
        u_i12, sigma_i12, lambda_i12);

[w3, w4, ...
    u_ic22, sigma_ic22, lambda_ic22, ...
    u_i12ic, sigma_i12ic, lambda_i12ic] ...
    = EM1_2_weighted(s2, Rs1ic, [w34/2, w34/2], ...
        u_ic2c, sigma_ic2c, lambda_ic2c, ...
        u_i12ic, sigma_i12ic, lambda_i12ic);


[w6, w7, ...
    u_ic2i1, sigma_ic2i1, lambda_ic2i1, ...
    u_i22, sigma_i22, lambda_i22] ...
    = EM1_2_weighted(s2, Rs1i1, [w567/2, w567/2], ...
        u_ic2c, sigma_ic2c, lambda_ic2c, ...
        u_i22, sigma_i22, lambda_i22);

w5 = 1e-4;
w7 = w7 - w5;

u_c2i1 = u_ic2i1 * 1.1;
sigma_c2i1 = sigma_ic2i1;
lambda_c2i1 = lambda_ic2i1;

ws = [w1, w2, w3, w4, w5, w6, w7];

Rs = rs_skewnorm(S, ws, ...
    [u_c1,        u_c1,        u_ic1,       u_ic1,        u_i11,       u_i11,        u_i11;
     u_ic2c,      u_i12c,      u_ic22,      u_i12ic,      u_c2i1,      u_ic2i1,      u_i22], ...
    [sigma_c1,    sigma_c1,    sigma_ic1,   sigma_ic1,    sigma_i11,   sigma_i11,    sigma_i11;
     sigma_ic2c,  sigma_i12c,  sigma_ic22,  sigma_i12ic,  sigma_c2i1,  sigma_ic2i1,  sigma_i22], ...
    [lambda_c1,   lambda_c1,   lambda_ic1,  lambda_ic1,   lambda_i11,  lambda_i11,   lambda_i11;
     lambda_ic2c, lambda_i12c, lambda_ic22, lambda_i12ic, lambda_c2i1, lambda_ic2i1, lambda_i22]);

ll = func_ll_skewnorm_mixture(S, ws, ...
    [u_c1,        u_c1,        u_ic1,       u_ic1,        u_i11,       u_i11,        u_i11;
     u_ic2c,      u_i12c,      u_ic22,      u_i12ic,      u_c2i1,      u_ic2i1,      u_i22], ...
    [sigma_c1,    sigma_c1,    sigma_ic1,   sigma_ic1,    sigma_i11,   sigma_i11,    sigma_i11;
     sigma_ic2c,  sigma_i12c,  sigma_ic22,  sigma_i12ic,  sigma_c2i1,  sigma_ic2i1,  sigma_i22], ...
    [lambda_c1,   lambda_c1,   lambda_ic1,  lambda_ic1,   lambda_i11,  lambda_i11,   lambda_i11;
     lambda_ic2c, lambda_i12c, lambda_ic22, lambda_i12ic, lambda_c2i1, lambda_ic2i1, lambda_i22]);

dll = ll - prev_ll;

i = 1;

fig_fit = figure('Position', [10,10,2000,1000]);
% ax_fit1 = subplot(2, 1, 1);
% ax_fit2 = subplot(2, 1, 2);

while dll > tolerance
    prev_ll = ll;
    
    plot_dist_10_xl_fn(fig_fit, S, '_2s_10j', ws, ...
        u_c1, sigma_c1, lambda_c1, ...
        u_ic1, sigma_ic1, lambda_ic1, ...
        u_i11, sigma_i11, lambda_i11, ...
        u_ic22, sigma_ic22, lambda_ic22, ...
        u_c2i1, sigma_c2i1, lambda_c2i1, ...
        u_ic2c, sigma_ic2c, lambda_ic2c, ...
        u_ic2i1, sigma_ic2i1, lambda_ic2i1, ...
        u_i12c, sigma_i12c, lambda_i12c, ...
        u_i12ic, sigma_i12ic, lambda_i12ic, ...
        u_i22, sigma_i22, lambda_i22)
    pause(0.001);
    
    Rs_old = Rs;
    Rs = rs_skewnorm(S, ws, ...
        [u_c1,        u_c1,        u_ic1,       u_ic1,        u_i11,       u_i11,        u_i11;
         u_ic2c,      u_i12c,      u_ic22,      u_i12ic,      u_c2i1,      u_ic2i1,      u_i22], ...
        [sigma_c1,    sigma_c1,    sigma_ic1,   sigma_ic1,    sigma_i11,   sigma_i11,    sigma_i11;
         sigma_ic2c,  sigma_i12c,  sigma_ic22,  sigma_i12ic,  sigma_c2i1,  sigma_ic2i1,  sigma_i22], ...
        [lambda_c1,   lambda_c1,   lambda_ic1,  lambda_ic1,   lambda_i11,  lambda_i11,   lambda_i11;
         lambda_ic2c, lambda_i12c, lambda_ic22, lambda_i12ic, lambda_c2i1, lambda_ic2i1, lambda_i22]);

    Rsj1 = Rs(1,:);
    Rsj2 = Rs(2,:);
    Rsj3 = Rs(3,:);
    Rsj4 = Rs(4,:);
    Rsj5 = Rs(5,:);
    Rsj6 = Rs(6,:);
    Rsj7 = Rs(7,:);

	w1 = sum(Rsj1) / M1;
	w2 = sum(Rsj2) / M1;
	w3 = sum(Rsj3) / M1;
	w4 = sum(Rsj4) / M1;
	w5 = sum(Rsj5) / M1;
	w6 = sum(Rsj6) / M1;
	w7 = sum(Rsj7) / M1;
    
    ws = [w1, w2, w3, w4, w5, w6, w7];
    
    [u_c1, sigma_c1, lambda_c1] = fit_skewnorm_weighted_one_step(s1, Rsj1+Rsj2, u_c1, sigma_c1, lambda_c1);
    [u_ic1, sigma_ic1, lambda_ic1] = fit_skewnorm_weighted_one_step(s1, Rsj3+Rsj4, u_ic1, sigma_ic1, lambda_ic1);
    [u_i11, sigma_i11, lambda_i11] = fit_skewnorm_weighted_one_step(s1, Rsj5+Rsj6+Rsj7, u_i11, sigma_i11, lambda_i11);
    [u_ic2c, sigma_ic2c, lambda_ic2c] = fit_skewnorm_weighted_one_step(s2, Rsj1, u_ic2c, sigma_ic2c, lambda_ic2c);
    [u_i12c, sigma_i12c, lambda_i12c] = fit_skewnorm_weighted_one_step(s2, Rsj2, u_i12c, sigma_i12c, lambda_i12c);
    [u_ic22, sigma_ic22, lambda_ic22] = fit_skewnorm_weighted_one_step(s2, Rsj3, u_ic22, sigma_ic22, lambda_ic22);
    [u_i12ic, sigma_i12ic, lambda_i12ic] = fit_skewnorm_weighted_one_step(s2, Rsj4, u_i12ic, sigma_i12ic, lambda_i12ic);
    [u_c2i1, sigma_c2i1, lambda_c2i1] = fit_skewnorm_weighted_one_step(s2, Rsj5, u_c2i1, sigma_c2i1, lambda_c2i1);
    [u_ic2i1, sigma_ic2i1, lambda_ic2i1] = fit_skewnorm_weighted_one_step(s2, Rsj6, u_ic2i1, sigma_ic2i1, lambda_ic2i1);
    [u_i22, sigma_i22, lambda_i22] = fit_skewnorm_weighted_one_step(s2, Rsj7, u_i22, sigma_i22, lambda_i22);
    
    ll = func_ll_skewnorm_mixture(S, ws, ...
        [u_c1,        u_c1,        u_ic1,       u_ic1,        u_i11,       u_i11,        u_i11;
         u_ic2c,      u_i12c,      u_ic22,      u_i12ic,      u_c2i1,      u_ic2i1,      u_i22], ...
        [sigma_c1,    sigma_c1,    sigma_ic1,   sigma_ic1,    sigma_i11,   sigma_i11,    sigma_i11;
         sigma_ic2c,  sigma_i12c,  sigma_ic22,  sigma_i12ic,  sigma_c2i1,  sigma_ic2i1,  sigma_i22], ...
        [lambda_c1,   lambda_c1,   lambda_ic1,  lambda_ic1,   lambda_i11,  lambda_i11,   lambda_i11;
         lambda_ic2c, lambda_i12c, lambda_ic22, lambda_i12ic, lambda_c2i1, lambda_ic2i1, lambda_i22]);
     
    dll = ll - prev_ll
    dlls(i) = dll;
    i = i + 1;
end


end
