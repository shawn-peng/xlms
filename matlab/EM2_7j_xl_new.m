function [dlls, alpha, beta, u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, u_i22, sigma_i22, lambda_i22] = EM2_4_xl(S,sl1,sl2,sl3,sl4)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
%   f2 = alpha*fi1 + beta*fc + (1-alpha-beta)*fi2
tollerance = 1e-8;
% tollerance = 1e-3;
S = S(:, S(2,:)~=0);
q1 = quantile(S(1,:), 0.01);
S = S(:, (S(1,:) > q1));
N = size(S, 1);

ll = 0;
prev_ll = -1;

s1 = S(1,:);
s2 = S(2,:);

M1 = size(s1, 2);
M2 = size(s2, 2);
 
S1_sorted = sort(s1, 'descend');
S2_sorted = sort(s2, 'descend');

[alpha, beta, ...
    u_c1, sigma_c1, lambda_c1, ...
    u_ic1, sigma_ic1, lambda_ic1, ...
    u_i11, sigma_i11, lambda_i11] = EM1_3i_xl(S, sl1, sl2, sl3);

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

% s2
u_c2 = u_c1;
sigma_c2 = sigma_c1;

[u_ic2, sigma_ic2] = fit_normal_weighted(s2, Rs1c);
[u_i12, sigma_i12] = fit_normal_weighted(s2, Rs1ic);
[u_i22, sigma_i22] = fit_normal_weighted(s2, Rs1i1);

lambda_c2  = sl1;
lambda_ic2 = sl1;
lambda_i12 = sl3;
lambda_i22 = sl4;

[u_ic2, sigma_ic2, lambda_ic2] = fit_skewnorm_weighted(s2, Rs1c,  u_ic2, sigma_ic2, lambda_ic2);
[u_i12, sigma_i12, lambda_i12] = fit_skewnorm_weighted(s2, Rs1ic, u_i12, sigma_i12, lambda_i12);
[u_i22, sigma_i22, lambda_i22] = fit_skewnorm_weighted(s2, Rs1i1, u_i22, sigma_i22, lambda_i22);


SIM_SIZE = 1000;

sim_xic2 = randn_skew([SIM_SIZE, 1], u_ic2, sigma_ic2, lambda_ic2);
sim_xi12 = randn_skew([SIM_SIZE, 1], u_i12, sigma_i12, lambda_i12);
sim_xi22 = randn_skew([SIM_SIZE, 1], u_i22, sigma_i22, lambda_i22);
p_ic2_i12 = sum(sim_xic2 > sim_xi12) / SIM_SIZE;
p_ic2_i22 = sum(sim_xic2 > sim_xi22) / SIM_SIZE;


w2 = 1e-24;
w1 = w1c - w2;
w3 = 1e-24;
w4 = 1e-24;
w5 = w1ic - w3;
w6 = 1e-24;
w7 = w1i1 - w4 - w6;

ws = [w1, w2, w3, w4, w5, w6, w7];


Rs = rs_skewnorm(S, ws, ...
    [u_c1,  u_c1,  u_ic1, u_i11, u_ic1, u_i11, u_i11;
     u_ic2, u_i12, u_c2,  u_c2,  u_i12, u_ic2, u_i22], ...
    [sigma_c1,  sigma_c1,  sigma_ic1, sigma_i11, sigma_ic1, sigma_i11, sigma_i11;
     sigma_ic2, sigma_i12, sigma_c2,  sigma_c2,  sigma_i12, sigma_ic2, sigma_i22], ...
    [lambda_c1,  lambda_c1,  lambda_ic1, lambda_i11, lambda_ic1, lambda_i11, lambda_i11;
     lambda_ic2, lambda_i12, lambda_c2,  lambda_c2,  lambda_i12, lambda_ic2, lambda_i22]);

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

dRs1 = Rsj1 - Rs1c;
dRs1(dRs1 < 0) = 0;

w_13 = w1 + w3;
w3 = sum(dRs1) / M1;
w1 = w_13 - w3;

[u_c2, sigma_c2, lambda_c2] = fit_skewnorm_weighted(s2, dRs1, u_c2, sigma_c2, -1);

rRs1 = Rsj1 - dRs1;
[w_12_1, w_12_2, ...
    u_ic2, sigma_ic2, lambda_ic2, ...
    u_i12c, sigma_i12c, lambda_i12c] ...
    = EM1_2_weighted(s2, rRs1, [w1/2, w1/2], ...
        u_ic2, sigma_ic2, lambda_ic2, ...
        u_i12, sigma_i12, lambda_i12);
% [u_ic2, sigma_ic2, lambda_ic2] = fit_skewnorm_weighted(s2, rRs1, u_ic2, sigma_ic2, 1);


ws = [w1, w2, w3, w4, w5, w6, w7];

fig = figure('Position', [10,10,2000,500]);
plot_dist_7_xl_fn(fig, S, '_s2_c7_xl', ws, ...
	u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, ...
	u_c2, sigma_c2, lambda_c2, u_ic2, sigma_ic2, lambda_ic2, u_i12, sigma_i12, lambda_i12, u_i22, sigma_i22, lambda_i22);
pause(.0005);

Rs = rs_skewnorm(S, ws, ...
    [u_c1,  u_c1,   u_ic1, u_i11, u_ic1, u_i11,   u_i11;
     u_ic2, u_i12c, u_c2,  u_c2,  u_i12, u_ic2, u_i22], ...
    [sigma_c1,  sigma_c1,   sigma_ic1, sigma_i11, sigma_ic1, sigma_i11,   sigma_i11;
     sigma_ic2, sigma_i12c, sigma_c2,  sigma_c2,  sigma_i12, sigma_ic2, sigma_i22], ...
    [lambda_c1,  lambda_c1,   lambda_ic1, lambda_i11, lambda_ic1, lambda_i11,   lambda_i11;
     lambda_ic2, lambda_i12c, lambda_c2,  lambda_c2,  lambda_i12, lambda_ic2, lambda_i22]);


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

% [w_125_1, w_125_2, ...
%     u_ic2, sigma_ic2, lambda_ic2, ...
%     u_i12, sigma_i12, lambda_i12] ...
%     = EM1_2_weighted(s2, Rsj1+Rsj2+Rsj5, [w1, w2+w5], ...
%         u_ic2, sigma_ic2, lambda_ic2, ...
%         u_i12, sigma_i12, lambda_i12);
% 
% w_12 = w1 + w2;
% w1 = w_125_1;
% w2 = w_12 - w_125_1;
% if w2 < 0
%     w2 = 1e-24;
% end
% w5 = w_125_2 - w2;
% assert(w2>0);

[w_67_1, w_67_2, ...
    u_ic2i1, sigma_ic2i1, lambda_ic2i1, ...
    u_i22, sigma_i22, lambda_i22] ...
    = EM1_2_weighted(s2, Rsj6+Rsj7, [(w6+w7)/2, (w6+w7)/2], ...
        u_ic2, sigma_ic2, lambda_ic2, ...
        u_i22, sigma_i22, lambda_i22);

w6 = w_67_1;
w7 = w_67_2;

ws = [w1, w2, w3, w4, w5, w6, w7];


fig = figure('Position', [10,10,2000,500]);
plot_dist_7_xl_fn(fig, S, '_s2_c7_xl', ws, ...
	u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, ...
	u_c2, sigma_c2, lambda_c2, u_ic2, sigma_ic2, lambda_ic2, u_i12, sigma_i12, lambda_i12, u_i22, sigma_i22, lambda_i22);
pause(.0005);

sc_fig = figure;
sc_ax = axes(sc_fig);

i = 1;
fig = figure('Position', [10,10,2000,500]);
while abs(ll - prev_ll) > tollerance
    dlls(i) = ll - prev_ll;
    i = i + 1;
    prev_ll = ll;

	% s1
    u_c1_new  = u_c1;
    u_ic1_new = u_ic1;
    u_i11_new = u_i11;
	% s2
    u_c2_new  = u_c2;
    u_ic2_new = u_ic2;
    u_i12_new = u_i12;
    u_i22_new = u_i22;

    sigma_c1_new  = sigma_c1;
    sigma_ic1_new = sigma_ic1;
    sigma_i11_new = sigma_i11;

    sigma_c2_new  = sigma_c2;
    sigma_ic2_new = sigma_ic2;
    sigma_i12_new = sigma_i12;
    sigma_i22_new = sigma_i22;

    lambda_c1_new  = lambda_c1;
    lambda_ic1_new = lambda_ic1;
    lambda_i11_new = lambda_i11;

    lambda_c2_new  = lambda_c2;
    lambda_ic2_new = lambda_ic2;
    lambda_i12_new = lambda_i12;
    lambda_i22_new = lambda_i22;

%     figure(fig);
	plot_dist_7_xl_fn(fig, S, '_s2_c7_xl', ws, ...
		u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, ...
		u_c2, sigma_c2, lambda_c2, u_ic2, sigma_ic2, lambda_ic2, u_i12, sigma_i12, lambda_i12, u_i22, sigma_i22, lambda_i22);
    pause(.0005);

    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
%     disp(ll);
%     disp(ll - prev_ll);
%     disp([alpha, beta, u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, u_i22, sigma_i22, lambda_i22]);
    
    
%     Rs_old = rsn_xl_d2i(s1, s2, ...
%         ws, ...
%         [u_c1,  u_c1,  u_ic1, u_i11, u_ic1, u_i11, u_i11;
%          u_ic2, u_i12, u_c2,  u_c2,  u_i12, u_ic2, u_i22], ...
%         [sigma_c1,  sigma_c1,  sigma_ic1, sigma_i11, sigma_ic1, sigma_i11, sigma_i11;
%          sigma_ic2, sigma_i12, sigma_c2,  sigma_c2,  sigma_i12, sigma_ic2, sigma_i22], ...
%         [lambda_c1,  lambda_c1,  lambda_ic1, lambda_i11, lambda_ic1, lambda_i11, lambda_i11;
%          lambda_ic2, lambda_i12, lambda_c2,  lambda_c2,  lambda_i12, lambda_ic2, lambda_i22]);

    Rs = rs_skewnorm(S, ws, ...
        [u_c1,  u_c1,   u_ic1, u_i11, u_ic1, u_i11,   u_i11;
         u_ic2, u_i12c, u_c2,  u_c2,  u_i12, u_ic2i1, u_i22], ...
        [sigma_c1,  sigma_c1,   sigma_ic1, sigma_i11, sigma_ic1, sigma_i11,   sigma_i11;
         sigma_ic2, sigma_i12c, sigma_c2,  sigma_c2,  sigma_i12, sigma_ic2i1, sigma_i22], ...
        [lambda_c1,  lambda_c1,   lambda_ic1, lambda_i11, lambda_ic1, lambda_i11,   lambda_i11;
         lambda_ic2, lambda_i12c, lambda_c2,  lambda_c2,  lambda_i12, lambda_ic2i1, lambda_i22]);

    Rsj1 = Rs(1,:);
    Rsj2 = Rs(2,:);
    Rsj3 = Rs(3,:);
    Rsj4 = Rs(4,:);
    Rsj5 = Rs(5,:);
    Rsj6 = Rs(6,:);
    Rsj7 = Rs(7,:);
    
    scatter(sc_ax, s1, s2, 10, Rsj1, '.');
    pause(0.001);

    Rs1c  = Rsj1 + Rsj2;
    Rs1ic = Rsj3 + Rsj5;
    Rs1i1 = Rsj4 + Rsj6 + Rsj7;

% 	Rs2c  = Rsj3 + Rsj4;
	Rs2c  = Rsj3;
	Rs2ic = Rsj1;
    Rs2ici1 = Rsj6;
	Rs2i1 = Rsj5;
    Rs2i1c = Rsj2;
	Rs2i2 = Rsj7;
    
    
    sum_Rs1c  = sum(Rs1c);
    sum_Rs1ic = sum(Rs1ic);
    sum_Rs1i1 = sum(Rs1i1);
    
    sum_Rs2c  = sum(Rs2c);
    sum_Rs2ic = sum(Rs2ic);
    sum_Rs2i1 = sum(Rs2i1);
    sum_Rs2i2 = sum(Rs2i2);

	w1 = sum(Rsj1) / M1;
	w2 = sum(Rsj2) / M1;
	w3 = sum(Rsj3) / M1;
	w4 = sum(Rsj4) / M1;
	w5 = sum(Rsj5) / M1;
	w6 = sum(Rsj6) / M1;
	w7 = sum(Rsj7) / M1;
    
    alpha_new = (sum_Rs1c + sum_Rs2c) / M1;
    beta_new = (sum_Rs1ic + sum_Rs2ic) / M1;
    
    ws_old = ws;
    ws = [w1, w2, w3, w4, w5, w6, w7];
    
    ll_old = func_ll2_7_xl(s1, s2, ws_old, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_c1_new, sigma_c1_new, lambda_c1_new] = fit_skewnorm_weighted_one_step(s1, Rs1c, u_c1, sigma_c1, lambda_c1);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_ic1_new, sigma_ic1_new, lambda_ic1_new] = fit_skewnorm_weighted_one_step(s1, Rs1ic, u_ic1, sigma_ic1, lambda_ic1);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_i11_new, sigma_i11_new, lambda_i11_new] = fit_skewnorm_weighted_one_step(s1, Rs1i1, u_i11, sigma_i11, lambda_i11);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_c2_new, sigma_c2_new, lambda_c2_new] = fit_skewnorm_weighted_one_step(s2, Rs2c, u_c2, sigma_c2, lambda_c2);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_ic2_new, sigma_ic2_new, lambda_ic2_new] = fit_skewnorm_weighted_one_step(s2, Rs2ic, u_ic2, sigma_ic2, lambda_ic2);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)

    [u_ic2i1_new, sigma_ic2i1_new, lambda_ic2i1_new] = fit_skewnorm_weighted_one_step(s2, Rs2ici1, u_ic2i1, sigma_ic2i1, lambda_ic2i1);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)

	[u_i12_new, sigma_i12_new, lambda_i12_new] = fit_skewnorm_weighted_one_step(s2, Rs2i1, u_i12, sigma_i12, lambda_i12);

    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_i12c_new, sigma_i12c_new, lambda_i12c_new] = fit_skewnorm_weighted_one_step(s2, Rs2i1c, u_i12c, sigma_i12c, lambda_i12c);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)
    
	[u_i22_new, sigma_i22_new, lambda_i22_new] = fit_skewnorm_weighted_one_step(s2, Rs2i2, u_i22, sigma_i22, lambda_i22);
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)

	u_c1  = u_c1_new;
	u_ic1 = u_ic1_new;
	u_i11 = u_i11_new;
	u_c2  = u_c2_new;
	u_ic2 = u_ic2_new;
	u_ic2i1 = u_ic2i1_new;
	u_i12 = u_i12_new;
	u_i12c = u_i12c_new;
	u_i22 = u_i22_new;

	sigma_c1  = sigma_c1_new;
	sigma_ic1 = sigma_ic1_new;
	sigma_i11 = sigma_i11_new;
	sigma_c2  = sigma_c2_new;
	sigma_ic2 = sigma_ic2_new;
	sigma_ic2i1 = sigma_ic2i1_new;
	sigma_i12 = sigma_i12_new;
	sigma_i12c = sigma_i12c_new;
	sigma_i22 = sigma_i22_new;

	lambda_c1  = lambda_c1_new;
	lambda_ic1 = lambda_ic1_new;
	lambda_i11 = lambda_i11_new;
	lambda_c2  = lambda_c2_new;
	lambda_ic2 = lambda_ic2_new;
	lambda_ic2i1 = lambda_ic2i1_new;
	lambda_i12 = lambda_i12_new;
	lambda_i12c = lambda_i12c_new;
	lambda_i22 = lambda_i22_new;
    
end

dlls = dlls(3:i-1);
end



