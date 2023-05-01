function [alpha, beta, u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, u_i22, sigma_i22, lambda_i22] = EM2_4_xl(S,sl1,sl2,sl3,sl4)
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
% s1 = s1(s1~=0);
% s2 = s2(s2~=0);

M1 = size(s1, 2);
M2 = size(s2, 2);
 
% S1_sorted = sort(S(1,:), 'descend');
% S2_sorted = sort(S(2,:), 'descend');
S1_sorted = sort(s1, 'descend');
S2_sorted = sort(s2, 'descend');
% S_sorted = S;

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

% % prior_thres = 38.5;
% prior_thres = S1_sorted(1) - (S1_sorted(1) - S1_sorted(M1)) * 1/3;
% prior_thres2 = S1_sorted(1) - (S1_sorted(1) - S1_sorted(M1)) * 3/4;
% nc1 = sum(S1_sorted > prior_thres);
% nic = sum(S1_sorted > prior_thres2 & S1_sorted <= prior_thres);
% nc2 = sum(S2_sorted > prior_thres);
% 
% alpha = nc1 / M1;
% % beta = (1-alpha) * 2/3; % assume the number of II is 2 times of IC
% beta = nic / M1;

% [u_c1, sigma_c1, lambda_c1] = sn_para_est(S1_sorted(1:nc1));
% % [u_ic1, sigma_ic1, lambda_ic1] = sn_para_est(S1_sorted(nc1:nic));
% [u_i11, sigma_i11, lambda_i11] = sn_para_est(S1_sorted(nic:M1));
% [u_i22, sigma_i22, lambda_i22] = sn_para_est(S2_sorted(int32(M2*(alpha+beta)):M2));

% [ndist, theta_c ] = fitdist_ntrunc(S1_sorted(1:nc1));
% [ndist, theta_ic] = fitdist_ntrunc(S1_sorted(nc1:nic));
% [ndist, theta_i1] = fitdist_ntrunc(S1_sorted(nic:M1));
% [ndist, theta_i2] = fitdist_ntrunc(S2_sorted(int32(M2*(alpha+beta)):M2));
% 
% [u_c1 , sigma_c1 ] = unpack2(theta_c );
% [u_ic1, sigma_ic1] = unpack2(theta_ic);
% [u_i11, sigma_i11] = unpack2(theta_i1);
% [u_i22, sigma_i22] = unpack2(theta_i2);

% [u_c1, sigma_c1] = normfit(S1_sorted(1:nc1));
% [u_ic1, sigma_ic1] = normfit(S1_sorted(nc1:nic));
% [u_i11, sigma_i11] = normfit(S1_sorted(nic:M1));
% [u_i22, sigma_i22] = normfit(S2_sorted(int32(M2*(alpha+beta)):M2));

% s2
u_c2 = u_c1;
sigma_c2 = sigma_c1;

[u_ic2, sigma_ic2] = fit_normal_weighted(s2, Rs1c);
[u_i12, sigma_i12] = fit_normal_weighted(s2, Rs1ic);
[u_i22, sigma_i22] = fit_normal_weighted(s2, Rs1i1);
% u_c2  = u_c1;
% u_ic2 = u_ic1;
% u_i12 = u_i11;
% 
% sigma_c2  = sigma_c1;
% sigma_ic2 = sigma_ic1;
% sigma_i12 = sigma_i11;

% lambda_c1  = 1;
% lambda_ic1 = 1;
% lambda_i11 = 1;

% lambda_c1 = 80;

% u_ic1 = (u_c1 + u_i11) / 2;
% sigma_ic1 = sigma_i11;

% lambda_c1 = lambda_c1 * sl1;
% lambda_ic1 = lambda_ic1 * sl2;
% lambda_i11 = lambda_i11 * sl3;

lambda_c2  = 1;
lambda_ic2 = sl2;
lambda_i12 = sl3;
lambda_i22 = sl4;

SIM_SIZE = 1000;

% sim_xc  = randn_skew([SIM_SIZE, 1], u_c1,  sigma_c1,  lambda_c1);
% sim_xic = randn_skew([SIM_SIZE, 1], u_ic1, sigma_ic1, lambda_ic1);
% sim_xi1 = randn_skew([SIM_SIZE, 1], u_i11, sigma_i11, lambda_i11);
% p_c_i1  = sum(sim_xc  > sim_xi1) / SIM_SIZE;
% p_c_ic  = sum(sim_xc  > sim_xic) / SIM_SIZE;
% p_ic_i1 = sum(sim_xic > sim_xi1) / SIM_SIZE;

% w1 = alpha * beta * p_c_ic;
% w2 = alpha * (1-beta) * p_c_i1;
% w3 = alpha * beta * (1-p_c_ic);
% w4 = alpha * (1-beta) * (1-p_c_i1);
% w5 = (1-alpha) * beta * p_ic_i1;
% w6 = (1-alpha) * beta * (1-p_ic_i1);
% w7 = (1-alpha) * (1-beta);

sim_xic2 = randn_skew([SIM_SIZE, 1], u_ic2, sigma_ic2, lambda_ic2);
sim_xi12 = randn_skew([SIM_SIZE, 1], u_i12, sigma_i12, lambda_i12);
sim_xi22 = randn_skew([SIM_SIZE, 1], u_i22, sigma_i22, lambda_i22);
p_ic2_i12 = sum(sim_xic2 > sim_xi12) / SIM_SIZE;
p_ic2_i22 = sum(sim_xic2 > sim_xi22) / SIM_SIZE;

% w1 = w1c * p_ic2_i12;
% w2 = w1c * (1 - p_ic2_i12);
% w3 = 1e-4;
% w4 = 1e-4;
% w5 = w1ic - w3;
% w6 = w1i1 * p_ic2_i22;
% w7 = w1i1 - w4 - w6;

w2 = 1e-4;
w1 = w1c - w2;
w3 = 1e-4;
w4 = 1e-4;
w5 = w1ic - w3;
w6 = 1e-4;
w7 = w1i1 - w4 - w6;

ws = [w1, w2, w3, w4, w5, w6, w7];

fig = figure('Position', [10,10,2000,500]);
plot_dist_7_xl_fn(S, '_s2_c7_xl', ws, ...
	u_c1, sigma_c1, lambda_c1, u_ic1, sigma_ic1, lambda_ic1, u_i11, sigma_i11, lambda_i11, ...
	u_c2, sigma_c2, lambda_c2, u_ic2, sigma_ic2, lambda_ic2, u_i12, sigma_i12, lambda_i12, u_i22, sigma_i22, lambda_i22);
pause(.0005);

fig = figure('Position', [10,10,2000,500]);
while abs(ll - prev_ll) > tollerance
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
	plot_dist_7_xl_fn(S, '_s2_c7_xl', ws, ...
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
    
    delta_c1 = lambda_c1 / sqrt(1+lambda_c1^2);
    Delta_c1 = sigma_c1 * delta_c1;
    Gamma_c1 = sigma_c1^2 * (1-delta_c1^2);
    delta_ic1 = lambda_ic1 / sqrt(1+lambda_ic1^2);
    Delta_ic1 = sigma_ic1 * delta_ic1;
    Gamma_ic1 = sigma_ic1^2 * (1-delta_ic1^2);
    delta_i11 = lambda_i11 / sqrt(1+lambda_i11^2);
    Delta_i11 = sigma_i11 * delta_i11;
    Gamma_i11 = sigma_i11^2 * (1-delta_i11^2);
    
    delta_c2 = lambda_c2 / sqrt(1+lambda_c2^2);
    Delta_c2 = sigma_c2 * delta_c2;
    Gamma_c2 = sigma_c2^2 * (1-delta_c2^2);
    delta_ic2 = lambda_ic2 / sqrt(1+lambda_ic2^2);
    Delta_ic2 = sigma_ic2 * delta_ic2;
    Gamma_ic2 = sigma_ic2^2 * (1-delta_ic2^2);
    delta_i12 = lambda_i12 / sqrt(1+lambda_i12^2);
    Delta_i12 = sigma_i12 * delta_i12;
    Gamma_i12 = sigma_i12^2 * (1-delta_i12^2);
    delta_i22 = lambda_i22 / sqrt(1+lambda_i22^2);
    Delta_i22 = sigma_i22 * delta_i22;
    Gamma_i22 = sigma_i22^2 * (1-delta_i22^2);

    [Vc1 , Wc1 ] = trunc_norm_moments(delta_c1  / sigma_c1  * (s1-u_c1 ), sqrt(1-delta_c1 ^2));
    [Vc2 , Wc2 ] = trunc_norm_moments(delta_c2  / sigma_c2  * (s2-u_c2 ), sqrt(1-delta_c2 ^2));
    [Vic1, Wic1] = trunc_norm_moments(delta_ic1 / sigma_ic1 * (s1-u_ic1), sqrt(1-delta_ic1^2));
    [Vic2, Wic2] = trunc_norm_moments(delta_ic2 / sigma_ic2 * (s2-u_ic2), sqrt(1-delta_ic2^2));
    [Vi11, Wi11] = trunc_norm_moments(delta_i11 / sigma_i11 * (s1-u_i11), sqrt(1-delta_i11^2));
    [Vi12, Wi12] = trunc_norm_moments(delta_i12 / sigma_i12 * (s2-u_i12), sqrt(1-delta_i12^2));
    [Vi22, Wi22] = trunc_norm_moments(delta_i22 / sigma_i22 * (s2-u_i22), sqrt(1-delta_i22^2));

    % sim_xc  = randn_skew([SIM_SIZE, 1], u_c1_new,  sigma_c1_new,  lambda_c1_new);
    % sim_xic = randn_skew([SIM_SIZE, 1], u_ic1_new, sigma_ic1_new, lambda_ic1_new);
    % sim_xi1 = randn_skew([SIM_SIZE, 1], u_i11_new, sigma_i11_new, lambda_i11_new);
    % p_c_i1  = sum(sim_xc  > sim_xi1) / SIM_SIZE;
    % p_c_ic  = sum(sim_xc  > sim_xic) / SIM_SIZE;
    % p_ic_i1 = sum(sim_xic > sim_xi1) / SIM_SIZE;
    
    Rs = rsn_xl_d2i(s1, s2, ...
        ws, ...
        [u_c1,  u_c1,  u_ic1, u_i11, u_ic1, u_i11, u_i11;
         u_ic2, u_i12, u_c2,  u_c2,  u_i12, u_ic2, u_i22], ...
        [sigma_c1,  sigma_c1,  sigma_ic1, sigma_i11, sigma_ic1, sigma_i11, sigma_i11;
         sigma_ic2, sigma_i12, sigma_c2,  sigma_c2,  sigma_i12, sigma_ic2, sigma_i22], ...
        [lambda_c1,  lambda_c1,  lambda_ic1, lambda_i11, lambda_ic1, lambda_i11, lambda_i11;
         lambda_ic2, lambda_i12, lambda_c2,  lambda_c2,  lambda_i12, lambda_ic2, lambda_i22]);
    
%     Rs_cic  = Rs(1,:);
%     Rs_ci1  = Rs(2,:);
%     Rs_ici1 = Rs(3,:);
%     Rs_i1i2 = Rs(4,:);
    
%     Rsj1 = Rs_cic * p_c_ic;
%     Rsj2 = Rs_ci1 * p_c_i1;
%     Rsj3 = Rs_cic * (1-p_c_ic);
%     Rsj4 = Rs_ci1 * (1-p_c_i1);
%     Rsj5 = Rs_ici1 * p_ic_i1;
%     Rsj6 = Rs_ici1 * (1-p_ic_i1);
%     Rsj7 = Rs_i1i2;
    
    Rsj1 = Rs(1,:);
    Rsj2 = Rs(2,:); 
    Rsj3 = Rs(3,:); 
    Rsj4 = Rs(4,:); 
    Rsj5 = Rs(5,:); 
    Rsj6 = Rs(6,:); 
    Rsj7 = Rs(7,:); 
    
    
%     Rs1c  = Rs_cic * p_c_ic     + Rs_ci1  * p_c_i1;
%     Rs1ic = Rs_cic * (1-p_c_ic) + Rs_ici1 * p_ic_i1;
%     Rs1i1 = Rs_i1i2             + Rs_ici1 * (1-p_ic_i1);
%     
%     Rs2c  = Rs_cic * (1-p_c_ic) + Rs_ci1  * (1-p_c_i1);
%     Rs2ic = Rs_cic * p_c_ic     + Rs_ici1 * (1-p_ic_i1);
%     Rs2i1 = Rs_ci1 * p_c_i1     + Rs_ici1 * p_ic_i1;
%     Rs2i2 = Rs_i1i2;

    Rs1c  = Rsj1 + Rsj2;
    Rs1ic = Rsj3 + Rsj5;
    Rs1i1 = Rsj4 + Rsj6 + Rsj7;

	Rs2c  = Rsj3 + Rsj4;
	Rs2ic = Rsj1 + Rsj6;
	Rs2i1 = Rsj2 + Rsj5;
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

%     alpha_new = w1 + w2 + w3 + w4;
%     beta_new = w1 + w3 + w5 + w6;
    
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
    

    u_c1_new  = (sum( Rs1c  .* (s1 - Vc1  * Delta_c1 ) )) / (sum_Rs1c);
    u_ic1_new = (sum( Rs1ic .* (s1 - Vic1 * Delta_ic1) )) / (sum_Rs1ic);
    u_i11_new = (sum( Rs1i1 .* (s1 - Vi11 * Delta_i11) )) / (sum_Rs1i1);

    u_c2_new  = (sum( Rs2c  .* (s2 - Vc2  * Delta_c2 ) )) / (sum_Rs2c);
    u_ic2_new = (sum( Rs2ic .* (s2 - Vic2 * Delta_ic2) )) / (sum_Rs2ic);
    u_i12_new = (sum( Rs2i1 .* (s2 - Vi12 * Delta_i12) )) / (sum_Rs2i1);
    u_i22_new = (sum( Rs2i2 .* (s2 - Vi22 * Delta_i22) )) / (sum_Rs2i2);

    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)

    Delta_c1_new  = (sum( Rs1c  .* Vc1  .* (s1 - u_c1_new ) ) ...
                    ) / (sum(Rs1c  .* Wc1 ));
    Delta_ic1_new = (sum( Rs1ic .* Vic1 .* (s1 - u_ic1_new) ) ...
                    ) / (sum(Rs1ic .* Wic1));
    Delta_i11_new = (sum( Rs1i1 .* Vi11 .* (s1 - u_i11_new) ) ...
                    ) / (sum(Rs1i1 .* Wi11));

    Delta_c2_new  = (sum( Rs2c  .* Vc2  .* (s2 - u_c2_new ) ) ...
                    ) / (sum(Rs2c  .* Wc2 ));
    Delta_ic2_new = (sum( Rs2ic .* Vic2 .* (s2 - u_ic2_new) ) ...
                    ) / (sum(Rs2ic .* Wic2 ));
    Delta_i12_new = (sum( Rs2i1 .* Vi12 .* (s2 - u_i12_new) ) ...
                    ) / (sum(Rs2i1 .* Wi12));
    Delta_i22_new = (sum( Rs2i2 .* Vi22 .* (s2 - u_i22_new) ) ...
                    ) / (sum(Rs2i2 .* Wi22));

    Gamma_c1_new  = ( sum( Rs1c .* ((s1 - u_c1_new).^2 - ...
                                  2 * Vc1 .* (s1 - u_c1_new) * Delta_c1_new + ...
                                  Wc1 * Delta_c1_new^2) ) ...
                     ) / (sum_Rs1c);
    Gamma_ic1_new = ( sum( Rs1ic .* ((s1 - u_ic1_new).^2 - ...
                                    2 * Vic1 .* (s1 - u_ic1_new) * Delta_ic1_new + ...
                                    Wic1 * Delta_ic1_new^2) ) ...
                     ) / (sum_Rs1ic);
    Gamma_i11_new = ( sum( Rs1i1 .* ((s1 - u_i11_new).^2 - ...
                                    2 * Vi11 .* (s1 - u_i11_new) * Delta_i11_new + ...
                                    Wi11 * Delta_i11_new^2) ) ...
                     ) / (sum_Rs1i1);

    Gamma_c2_new  = ( sum( Rs2c .* ((s2 - u_c2_new).^2 - ...
                                  2 * Vc2 .* (s2 - u_c2_new) * Delta_c2_new + ...
                                  Wc2 * Delta_c2_new^2) ) ...
                     ) / (sum_Rs2c);
    Gamma_ic2_new = ( sum( Rs2ic .* ((s2 - u_ic2_new).^2 - ...
                                    2 * Vic2 .* (s2 - u_ic2_new) * Delta_ic2_new + ...
                                    Wic2 * Delta_ic2_new^2) ) ...
                     ) / (sum_Rs2ic);
    Gamma_i12_new = ( sum( Rs2i1 .* ((s2 - u_i12_new).^2 - ...
                                    2 * Vi12 .* (s2 - u_i12_new) * Delta_i12_new + ...
                                    Wi12 * Delta_i12_new^2) ) ...
                     ) / (sum_Rs2i1);
    Gamma_i22_new = ( sum( Rs2i2 .* ((s2 - u_i22_new).^2 - ...
                                    2 * Vi22 .* (s2 - u_i22_new) * Delta_i22_new + ...
                                    Wi22 * Delta_i22_new^2) ) ...
                     ) / (sum_Rs2i2);

    Delta_c1 = Delta_c1_new;
    Delta_ic1 = Delta_ic1_new;
    Delta_i11 = Delta_i11_new;
    Delta_c2 = Delta_c2_new;
    Delta_ic2 = Delta_ic2_new;
    Delta_i12 = Delta_i12_new;
    Delta_i22 = Delta_i22_new;

    Gamma_c1 = Gamma_c1_new;
    Gamma_ic1 = Gamma_ic1_new;
    Gamma_i11 = Gamma_i11_new;
    Gamma_c2 = Gamma_c2_new;
    Gamma_ic2 = Gamma_ic2_new;
    Gamma_i12 = Gamma_i12_new;
    Gamma_i22 = Gamma_i22_new;

    lambda_c1 = sign(Delta_c1) * sqrt(Delta_c1^2 / Gamma_c1);
    lambda_ic1 = sign(Delta_ic1) * sqrt(Delta_ic1^2 / Gamma_ic1);
    lambda_i11 = sign(Delta_i11) * sqrt(Delta_i11^2 / Gamma_i11);
    lambda_c2 = sign(Delta_c2) * sqrt(Delta_c2^2 / Gamma_c2);
    lambda_ic2 = sign(Delta_ic2) * sqrt(Delta_ic2^2 / Gamma_ic2);
    lambda_i12 = sign(Delta_i12) * sqrt(Delta_i12^2 / Gamma_i12);
    lambda_i22 = sign(Delta_i22) * sqrt(Delta_i22^2 / Gamma_i22);

    lambda_c1_new  = lambda_c1;
    lambda_ic1_new  = lambda_ic1;
    lambda_i11_new = lambda_i11;
    lambda_c2_new  = lambda_c2;
    lambda_ic2_new  = lambda_ic2;
    lambda_i12_new = lambda_i12;
    lambda_i22_new = lambda_i22;
    
    sigma_c1 = sqrt(Gamma_c1 + Delta_c1^2);
    sigma_ic1 = sqrt(Gamma_ic1 + Delta_ic1^2);
    sigma_i11 = sqrt(Gamma_i11 + Delta_i11^2);
    sigma_c2 = sqrt(Gamma_c2 + Delta_c2^2);
    sigma_ic2 = sqrt(Gamma_ic2 + Delta_ic2^2);
    sigma_i12 = sqrt(Gamma_i12 + Delta_i12^2);
    sigma_i22 = sqrt(Gamma_i22 + Delta_i22^2);

    sigma_c1_new  = sigma_c1;
    sigma_ic1_new = sigma_ic1;
    sigma_i11_new = sigma_i11;
    sigma_c2_new  = sigma_c2;
    sigma_ic2_new = sigma_ic2;
    sigma_i12_new = sigma_i12;
    sigma_i22_new = sigma_i22;
    
    
    ll = func_ll2_7_xl(s1, s2, ws, ...
        u_c1_new, sigma_c1_new, lambda_c1_new, u_ic1_new, sigma_ic1_new, lambda_ic1_new, u_i11_new, sigma_i11_new, lambda_i11_new, ...
        u_c2_new, sigma_c2_new, lambda_c2_new, u_ic2_new, sigma_ic2_new, lambda_ic2_new, u_i12_new, sigma_i12_new, lambda_i12_new, ...
        u_i22_new, sigma_i22_new, lambda_i22_new)

    
%     lambda_1 = (2*M1 + sum_Rs1i1 + 2*M2 - sum_Rs2i2);
%     w1 = p_c_ic * (sum_Rs1c + sum_Rs2ic + sum_Rs1ic + sum_Rs2c) / lambda_1;
%     w2 = p_c_i1 * (sum_Rs1c + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2c) / lambda_1;
%     w3 = (1-p_c_ic) * (sum_Rs1c + sum_Rs2ic + sum_Rs1ic + sum_Rs2c) / lambda_1;
%     w4 = (1-p_c_i1) * (sum_Rs1c + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2c) / lambda_1;
%     w5 = p_ic_i1 * (sum_Rs1ic + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2ic) / lambda_1;
%     w6 = (1-p_ic_i1) * (sum_Rs1ic + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2ic) / lambda_1;
%     w7 = (sum_Rs1i1 + sum_Rs2i2) / lambda_1;


    
    alpha = alpha_new;
    beta = beta_new;
    u_c1 = u_c1_new;
    u_ic1 = u_ic1_new;
    u_i11 = u_i11_new;
    u_c2 = u_c2_new;
    u_ic2 = u_ic2_new;
    u_i12 = u_i12_new;
    u_i22 = u_i22_new;

    Delta_c1 = Delta_c1_new;
    Delta_ic1 = Delta_ic1_new;
    Delta_i11 = Delta_i11_new;
    Delta_c2 = Delta_c2_new;
    Delta_ic2 = Delta_ic2_new;
    Delta_i12 = Delta_i12_new;
    Delta_i22 = Delta_i22_new;

    Gamma_c1 = Gamma_c1_new;
    Gamma_ic1 = Gamma_ic1_new;
    Gamma_i11 = Gamma_i11_new;
    Gamma_c2 = Gamma_c2_new;
    Gamma_ic2 = Gamma_ic2_new;
    Gamma_i12 = Gamma_i12_new;
    Gamma_i22 = Gamma_i22_new;

    lambda_c1 = sign(Delta_c1) * sqrt(Delta_c1^2 / Gamma_c1);
    lambda_ic1 = sign(Delta_ic1) * sqrt(Delta_ic1^2 / Gamma_ic1);
    lambda_i11 = sign(Delta_i11) * sqrt(Delta_i11^2 / Gamma_i11);
    lambda_c2 = sign(Delta_c2) * sqrt(Delta_c2^2 / Gamma_c2);
    lambda_ic2 = sign(Delta_ic2) * sqrt(Delta_ic2^2 / Gamma_ic2);
    lambda_i12 = sign(Delta_i12) * sqrt(Delta_i12^2 / Gamma_i12);
    lambda_i22 = sign(Delta_i22) * sqrt(Delta_i22^2 / Gamma_i22);

    sigma_c1 = sqrt(Gamma_c1 + Delta_c1^2);
    sigma_ic1 = sqrt(Gamma_ic1 + Delta_ic1^2);
    sigma_i11 = sqrt(Gamma_i11 + Delta_i11^2);
    sigma_c2 = sqrt(Gamma_c2 + Delta_c2^2);
    sigma_ic2 = sqrt(Gamma_ic2 + Delta_ic2^2);
    sigma_i12 = sqrt(Gamma_i12 + Delta_i12^2);
    sigma_i22 = sqrt(Gamma_i22 + Delta_i22^2);
    
%     disp(ll - prev_ll);
%     disp(alpha);
%     sigma_c1
%     sigma_i
%     u_c1
%     u_i
%     lambda_c1
%     lambda_i
end

end



    
%     syms sym_alpha sym_beta;
%     syms sym_w1c(sym_alpha, sym_beta) sym_w1ic(sym_alpha, sym_beta) sym_wi1(sym_alpha, sym_beta);
%     syms sym_w2c(sym_alpha, sym_beta) sym_w2ic(sym_alpha, sym_beta) sym_w2i1(sym_alpha, sym_beta) sym_w2i2(sym_alpha, sym_beta);
%     
%     sym_w1c  = sym_alpha * sym_beta * p_c_ic     + sym_alpha     * (1-sym_beta) * p_c_i1;
%     sym_w1ic = sym_alpha * sym_beta * (1-p_c_ic) + (1-sym_alpha) * sym_beta     * p_ic_i1;
%     sym_w1i1 = 1 - sym_w1c - sym_w1ic;
% 
%     sym_w2c  = sym_alpha * sym_beta * (1-p_c_ic) + sym_alpha     * (1-sym_beta) * (1-p_c_i1);
%     sym_w2ic = sym_alpha * sym_beta * p_c_ic     + (1-sym_alpha) * sym_beta     * (1-p_ic_i1);
%     sym_w2i2 = (1-sym_alpha) * (1-sym_beta);
%     sym_w2i1 = 1 - sym_w2c - sym_w2ic - sym_w2i2;
    
%     w1_vec = [sym_w1c; sym_w1ic; sym_w1i1];
%     w2_vec = [sym_w2c; sym_w2ic; sym_w2i2; sym_w2i1];
    
%     d_w1_d_alpha = diff(w1_vec, sym_alpha);
%     d_w1_d_beta  = diff(w1_vec, sym_beta);
%     d_w2_d_alpha = diff(w2_vec, sym_alpha);
%     d_w2_d_beta  = diff(w2_vec, sym_beta);
%     
%     sum_Rs1 = [sum_Rs1c; sum_Rs1ic; sum_Rs1i1];
%     sum_Rs2 = [sum_Rs2c; sum_Rs2ic; sum_Rs2i2; sum_Rs2i1];
%     sum_Rs1 = sum(Rs1, 2);
%     sum_Rs2 = sum(Rs2, 2);
    
%     eqns = simplify([
%         d_w1_d_alpha' * (sum_Rs1 ./ w1_vec)
%             + d_w2_d_alpha' * (sum_Rs2 ./ w2_vec) == 0;
%         d_w1_d_beta'  * (sum_Rs1 ./ w1_vec)
%             + d_w2_d_beta'  * (sum_Rs2 ./ w2_vec) == 0;
%     ]);
%     
%     eqns = simplify([
%         d_w1_d_alpha' * (sum_Rs1 ./ w1_vec)
%             + d_w2_d_alpha' * (sum_Rs2 ./ w2_vec) == 0;
%         d_w1_d_beta'  * (sum_Rs1 ./ w1_vec)
%             + d_w2_d_beta'  * (sum_Rs2 ./ w2_vec) == 0;
%     ]);

%     assume([
%         sym_alpha >= 0;
%         sym_alpha <= 1;
%         sym_beta  >= 0;
%         sym_beta  <= 1;
%     ]);

%     E = solve(eqns, [sym_alpha, sym_beta], 'IgnoreAnalyticConstraints',1)
    
%     eqns = [
%         sym_w1c  == sym_alpha * sym_beta * p_c_ic     + sym_alpha     * (1-sym_beta) * p_c_i1;
%         sym_w1ic == sym_alpha * sym_beta * (1-p_c_ic) + (1-sym_alpha) * sym_beta     * p_ic_i1;
%         % sym_w1i1

%         sym_w2c  == sym_alpha * sym_beta * (1-p_c_ic) + sym_alpha     * (1-sym_beta) * (1-p_c_i1);
%         sym_w2ic == sym_alpha * sym_beta * p_c_ic     + (1-sym_alpha) * sym_beta     * (1-p_ic_i1);
%         % sym_w2i1
%         sym_w2i2 == (1-sym_alpha) * (1-sym_beta);
%     ];
%     E = solve(eqns, [sym_alpha sym_beta]);
