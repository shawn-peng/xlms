function [alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2] = EM2_4_xl(S,sl1,sl2,sl3,sl4)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
%   f2 = alpha*fi1 + beta*fc + (1-alpha-beta)*fi2
tollerance = 1e-8;
% tollerance = 1e-3;
N = size(S,1);

ll = 0;
prev_ll = -1;

s1 = S(1,:);
s1 = s1(s1~=0);
s2 = S(2,:);
s2 = s2(s2~=0);

M1 = size(s1,2);
M2 = size(s2,2);
 
% S1_sorted = sort(S(1,:), 'descend');
% S2_sorted = sort(S(2,:), 'descend');
S1_sorted = sort(s1, 'descend');
S2_sorted = sort(s2, 'descend');
% S_sorted = S;

% prior_thres = 38.5;
prior_thres = S1_sorted(1) - (S1_sorted(1) - S1_sorted(M1))/3;
prior_thres2 = S1_sorted(1) - (S1_sorted(1) - S1_sorted(M1)) * 2/3;
nc1 = sum(S1_sorted > prior_thres);
nic = sum(S1_sorted > prior_thres2 & S1_sorted <= prior_thres);
nc2 = sum(S2_sorted > prior_thres);

alpha = nc1 / M1;
% beta = (1-alpha) * 2/3; % assume the number of II is 2 times of IC
beta = nic / M1;

% [u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:nc1));
% % [u_ic, sigma_ic, lambda_ic] = sn_para_est(S1_sorted(nc1:nic));
% [u_i1, sigma_i1, lambda_i1] = sn_para_est(S1_sorted(nic:M1));
% [u_i2, sigma_i2, lambda_i2] = sn_para_est(S2_sorted(int32(M2*(alpha+beta)):M2));

% [ndist, theta_c ] = fitdist_ntrunc(S1_sorted(1:nc1));
% [ndist, theta_ic] = fitdist_ntrunc(S1_sorted(nc1:nic));
% [ndist, theta_i1] = fitdist_ntrunc(S1_sorted(nic:M1));
% [ndist, theta_i2] = fitdist_ntrunc(S2_sorted(int32(M2*(alpha+beta)):M2));
% 
% [u_c , sigma_c ] = unpack2(theta_c );
% [u_ic, sigma_ic] = unpack2(theta_ic);
% [u_i1, sigma_i1] = unpack2(theta_i1);
% [u_i2, sigma_i2] = unpack2(theta_i2);

[u_c, sigma_c] = normfit(S1_sorted(1:nc1));
[u_ic, sigma_ic] = normfit(S1_sorted(nc1:nic));
[u_i1, sigma_i1] = normfit(S1_sorted(nic:M1));
[u_i2, sigma_i2] = normfit(S2_sorted(int32(M2*(alpha+beta)):M2));

lambda_c  = 1;
lambda_ic = 1;
lambda_i1 = 1;
lambda_i2 = 1;

% lambda_c = 80;

% u_ic = (u_c + u_i1) / 2;
% sigma_ic = sigma_i1;

lambda_c = lambda_c * sl1;
lambda_ic = lambda_ic * sl2;
lambda_i1 = lambda_i1 * sl3;
lambda_i2 = lambda_i2 * sl4;

SIM_SIZE = 1000;

sim_xc  = randn_skew([SIM_SIZE, 1], u_c,  sigma_c,  lambda_c);
sim_xic = randn_skew([SIM_SIZE, 1], u_ic, sigma_ic, lambda_ic);
sim_xi1 = randn_skew([SIM_SIZE, 1], u_i1, sigma_i1, lambda_i1);
p_c_i1  = sum(sim_xc  > sim_xi1) / SIM_SIZE;
p_c_ic  = sum(sim_xc  > sim_xic) / SIM_SIZE;
p_ic_i1 = sum(sim_xic > sim_xi1) / SIM_SIZE;

w1 = alpha * beta * p_c_ic;
w2 = alpha * (1-beta) * p_c_i1;
w3 = alpha * beta * (1-p_c_ic);
w4 = alpha * (1-beta) * (1-p_c_i1);
w5 = (1-alpha) * beta * p_ic_i1;
w6 = (1-alpha) * beta * (1-p_ic_i1);
w7 = (1-alpha) * (1-beta);

ws = [w1, w2, w3, w4, w5, w6, w7];

fig = figure('Position', [10,10,2000,500]);
plot_dist_xl_fn(S, '_s2_c4_xl', ws, ...
    u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2);
pause(.0005);

fig = figure('Position', [10,10,2000,500]);
while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

    u_c_new  = u_c;
    u_ic_new = u_ic;
    u_i1_new = u_i1;
    u_i2_new = u_i2;

    sigma_c_new  = sigma_c;
    sigma_ic_new = sigma_ic;
    sigma_i1_new = sigma_i1;
    sigma_i2_new = sigma_i2;

    lambda_c_new  = lambda_c;
    lambda_ic_new = lambda_ic;
    lambda_i1_new = lambda_i1;
    lambda_i2_new = lambda_i2;

%     figure(fig);
    plot_dist_xl_fn(S, '_s2_c4_xl', ws, ...
        u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
        u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2);
    pause(.0005);

    ll = func_ll2_4_xl(s1, s2, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new)
%     disp(ll);
%     disp(ll - prev_ll);
%     disp([alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2]);
    
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);
    delta_ic = lambda_ic / sqrt(1+lambda_ic^2);
    Delta_ic = sigma_ic * delta_ic;
    Gamma_ic = sigma_ic^2 * (1-delta_ic^2);
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);
    delta_i2 = lambda_i2 / sqrt(1+lambda_i2^2);
    Delta_i2 = sigma_i2 * delta_i2;
    Gamma_i2 = sigma_i2^2 * (1-delta_i2^2);

    [Vc1 , Wc1 ] = trunc_norm_moments(delta_c  / sigma_c  * (s1-u_c ), sqrt(1-delta_c ^2));
    [Vc2 , Wc2 ] = trunc_norm_moments(delta_c  / sigma_c  * (s2-u_c ), sqrt(1-delta_c ^2));
    [Vic1, Wic1] = trunc_norm_moments(delta_ic / sigma_ic * (s1-u_ic), sqrt(1-delta_ic^2));
    [Vic2, Wic2] = trunc_norm_moments(delta_ic / sigma_ic * (s2-u_ic), sqrt(1-delta_ic^2));
    [Vi11, Wi11] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));
    [Vi12, Wi12] = trunc_norm_moments(delta_i1 / sigma_i1 * (s2-u_i1), sqrt(1-delta_i1^2));
    [Vi22, Wi22] = trunc_norm_moments(delta_i2 / sigma_i2 * (s2-u_i2), sqrt(1-delta_i2^2));
    
%     p_s1_c  = alpha * beta * p_c_ic     + alpha     * (1-beta) * p_c_i1;
%     p_s1_ic = alpha * beta * (1-p_c_ic) + (1-alpha) * beta     * p_ic_i1;
%     p_s1_i1 = 1 - p_s1_c - p_s1_ic;
%     p_s2_c  = alpha * beta * (1-p_c_ic) + alpha     * (1-beta) * (1-p_c_i1);
%     p_s2_ic = alpha * beta * p_c_ic     + (1-alpha) * beta     * (1-p_ic_i1);
%     p_s2_i2 = (1-alpha) * (1-beta);
%     p_s2_i1 = 1 - p_s2_c - p_s2_ic - p_s2_i2;

    p_s1_c = w1 + w2;
    p_s1_ic = w3 + w5;
    p_s1_i1 = w4 + w6 + w7;
    p_s2_c = w3 + w4;
    p_s2_ic = w1 + w6;
    p_s2_i1 = w2 + w5;
    p_s2_i2 = w7;
    
    Rs1 = rsn(s1, ...
        [p_s1_c, p_s1_ic, p_s1_i1], ...
        [u_c, u_ic, u_i1], [sigma_c, sigma_ic, sigma_i1], ...
        [lambda_c, lambda_ic, lambda_i1]);
    Rs2 = rsn(s2, ...
        [p_s2_c, p_s2_ic, p_s2_i1, p_s2_i2], ...
        [u_c, u_ic, u_i1, u_i2], ...
        [sigma_c, sigma_ic, sigma_i1, sigma_i2], ...
        [lambda_c, lambda_ic, lambda_i1, lambda_i2]);
    
    Rs1c  = Rs1(1,:);
    Rs1ic = Rs1(2,:);
    Rs1i1 = Rs1(3,:);
    
    Rs2c  = Rs2(1,:);
    Rs2ic = Rs2(2,:);
    Rs2i1 = Rs2(3,:);
    Rs2i2 = Rs2(4,:);
    
    sum_Rs1c  = sum(Rs1c);
    sum_Rs1ic = sum(Rs1ic);
    sum_Rs1i1 = sum(Rs1i1);
    
    sum_Rs2c  = sum(Rs2c);
    sum_Rs2ic = sum(Rs2ic);
    sum_Rs2i1 = sum(Rs2i1);
    sum_Rs2i2 = sum(Rs2i2);
    

    sim_xc  = randn_skew([SIM_SIZE, 1], u_c_new,  sigma_c_new,  lambda_c_new);
    sim_xic = randn_skew([SIM_SIZE, 1], u_ic_new, sigma_ic_new, lambda_ic_new);
    sim_xi1 = randn_skew([SIM_SIZE, 1], u_i1_new, sigma_i1_new, lambda_i1_new);
    p_c_i1  = sum(sim_xc  > sim_xi1) / SIM_SIZE;
    p_c_ic  = sum(sim_xc  > sim_xic) / SIM_SIZE;
    p_ic_i1 = sum(sim_xic > sim_xi1) / SIM_SIZE;
    
    lambda_1 = (2*M1 + sum_Rs1i1 + 2*M2 - sum_Rs2i2);
    w1 = p_c_ic * (sum_Rs1c + sum_Rs2ic + sum_Rs1ic + sum_Rs2c) / lambda_1;
    w2 = p_c_i1 * (sum_Rs1c + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2c) / lambda_1;
    w3 = (1-p_c_ic) * (sum_Rs1c + sum_Rs2ic + sum_Rs1ic + sum_Rs2c) / lambda_1;
    w4 = (1-p_c_i1) * (sum_Rs1c + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2c) / lambda_1;
    w5 = p_ic_i1 * (sum_Rs1ic + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2ic) / lambda_1;
    w6 = (1-p_ic_i1) * (sum_Rs1ic + sum_Rs2i1 + sum_Rs1i1 + sum_Rs2ic) / lambda_1;
    w7 = (sum_Rs1i1 + sum_Rs2i2) / lambda_1;
    
    ws = [w1, w2, w3, w4, w5, w6, w7];
    
    ll = func_ll2_4_xl(s1, s2, ws, ...
    u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
    u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new)



    u_c_new  = (sum( Rs1c  .* (s1 - Vc1  * Delta_c ) ) + ...
                sum( Rs2c  .* (s2 - Vc2  * Delta_c ) )) / (sum_Rs1c + sum_Rs2c);
    u_ic_new = (sum( Rs1ic .* (s1 - Vic1 * Delta_ic) ) + ...
                sum( Rs2ic .* (s2 - Vic2 * Delta_ic) )) / (sum_Rs1ic + sum_Rs2ic);
    u_i1_new = (sum( Rs1i1 .* (s1 - Vi11 * Delta_i1) ) + ...
                sum( Rs2i1 .* (s2 - Vi12 * Delta_i1) )) / (sum_Rs1i1 + sum_Rs2i1);
    u_i2_new = (sum( Rs2i2 .* (s2 - Vi22 * Delta_i2) )) / (sum_Rs2i2);

    ll = func_ll2_4_xl(s1, s2, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new)

    Delta_c_new  = (sum( Rs1c  .* Vc1  .* (s1 - u_c_new ) ) + ...
                    sum( Rs2c  .* Vc2  .* (s2 - u_c_new ) ) ...
                    ) / (sum(Rs1c  .* Wc1 ) + sum(Rs2c  .* Wc2 ));
    Delta_ic_new = (sum( Rs1ic .* Vic1 .* (s1 - u_ic_new) ) + ...
                    sum( Rs2ic .* Vic2 .* (s2 - u_ic_new) ) ...
                    ) / (sum(Rs1ic .* Wic1) + sum(Rs2ic .* Wic2 ));
    Delta_i1_new = (sum( Rs1i1 .* Vi11 .* (s1 - u_i1_new) ) + ...
                    sum( Rs2i1 .* Vi12 .* (s2 - u_i1_new) ) ...
                    ) / (sum(Rs1i1 .* Wi11) + sum(Rs2i1 .* Wi12));
    Delta_i2_new = (sum( Rs2i2 .* Vi22 .* (s2 - u_i2_new) ) ...
                    ) / (sum(Rs2i2 .* Wi22));

    Gamma_c_new  = ( sum( Rs1c .* ((s1 - u_c_new).^2 - ...
                                  2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + ...
                                  Wc1 * Delta_c_new^2) ) + ...
                     sum( Rs2c .* ((s2 - u_c_new).^2 - ...
                                  2 * Vc2 .* (s2 - u_c_new) * Delta_c_new + ...
                                  Wc2 * Delta_c_new^2) ) ...
                     ) / (sum_Rs1c + sum_Rs2c);
    Gamma_ic_new = ( sum( Rs1ic .* ((s1 - u_ic_new).^2 - ...
                                    2 * Vic1 .* (s1 - u_ic_new) * Delta_ic_new + ...
                                    Wic1 * Delta_ic_new^2) ) + ...
                     sum( Rs2ic .* ((s2 - u_ic_new).^2 - ...
                                    2 * Vic2 .* (s2 - u_ic_new) * Delta_ic_new + ...
                                    Wic2 * Delta_ic_new^2) ) ...
                     ) / (sum_Rs1ic + sum_Rs2ic);
    Gamma_i1_new = ( sum( Rs1i1 .* ((s1 - u_i1_new).^2 - ...
                                    2 * Vi11 .* (s1 - u_i1_new) * Delta_i1_new + ...
                                    Wi11 * Delta_i1_new^2) ) + ...
                     sum( Rs2i1 .* ((s2 - u_i1_new).^2 - ...
                                    2 * Vi12 .* (s2 - u_i1_new) * Delta_i1_new + ...
                                    Wi12 * Delta_i1_new^2) ) ...
                     ) / (sum_Rs1i1 + sum_Rs2i1);
    Gamma_i2_new = ( sum( Rs2i2 .* ((s2 - u_i2_new).^2 - ...
                                    2 * Vi22 .* (s2 - u_i2_new) * Delta_i2_new + ...
                                    Wi22 * Delta_i2_new^2) ) ...
                     ) / (sum_Rs2i2);

    Delta_c = Delta_c_new;
    Delta_ic = Delta_ic_new;
    Delta_i1 = Delta_i1_new;
    Delta_i2 = Delta_i2_new;
    Gamma_c = Gamma_c_new;
    Gamma_ic = Gamma_ic_new;
    Gamma_i1 = Gamma_i1_new;
    Gamma_i2 = Gamma_i2_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_ic = sign(Delta_ic) * sqrt(Delta_ic^2 / Gamma_ic);
    lambda_i1 = sign(Delta_i1) * sqrt(Delta_i1^2 / Gamma_i1);
    lambda_i2 = sign(Delta_i2) * sqrt(Delta_i2^2 / Gamma_i2);
    lambda_c_new  = lambda_c;
    lambda_ic_new  = lambda_ic;
    lambda_i1_new = lambda_i1;
    lambda_i2_new = lambda_i2;
    
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    sigma_ic = sqrt(Gamma_ic + Delta_ic^2);
    sigma_i1 = sqrt(Gamma_i1 + Delta_i1^2);
    sigma_i2 = sqrt(Gamma_i2 + Delta_i2^2);
    sigma_c_new  = sigma_c;
    sigma_ic_new = sigma_ic;
    sigma_i1_new = sigma_i1;
    sigma_i2_new = sigma_i2;
    
    
    ll = func_ll2_4_xl(s1, s2, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new)



    alpha_new = w1 + w2 + w3 + w4;
    beta_new = w1 + w3 + w5 + w6;

    
    alpha = alpha_new;
    beta = beta_new;
    u_c = u_c_new;
    u_ic = u_ic_new;
    u_i1 = u_i1_new;
    u_i2 = u_i2_new;
    Delta_c = Delta_c_new;
    Delta_ic = Delta_ic_new;
    Delta_i1 = Delta_i1_new;
    Delta_i2 = Delta_i2_new;
    Gamma_c = Gamma_c_new;
    Gamma_ic = Gamma_ic_new;
    Gamma_i1 = Gamma_i1_new;
    Gamma_i2 = Gamma_i2_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_ic = sign(Delta_ic) * sqrt(Delta_ic^2 / Gamma_ic);
    lambda_i1 = sign(Delta_i1) * sqrt(Delta_i1^2 / Gamma_i1);
    lambda_i2 = sign(Delta_i2) * sqrt(Delta_i2^2 / Gamma_i2);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    sigma_ic = sqrt(Gamma_ic + Delta_ic^2);
    sigma_i1 = sqrt(Gamma_i1 + Delta_i1^2);
    sigma_i2 = sqrt(Gamma_i2 + Delta_i2^2);
    
%     disp(ll - prev_ll);
%     disp(alpha);
%     sigma_c
%     sigma_i
%     u_c
%     u_i
%     lambda_c
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
