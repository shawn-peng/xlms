function [alpha, beta, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1] = EM1_3ic_xl(S,sl1,sl2,sl3,varargin)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
%   f2 = alpha*fi1 + beta*fc + (1-alpha-beta)*fi2

parser = inputParser;
addParameter(parser, 'tolerance', 1e-7);
addParameter(parser, 'lower_quantile', 0.01);
addParameter(parser, 'c_range', 1/2);
addParameter(parser, 'constraints', []);

parse(parser, varargin{:});
tolerance       = parser.Results.tolerance;
lower_quantile  = parser.Results.lower_quantile;
c_range         = parser.Results.c_range;
constraints     = parser.Results.constraints;

q1 = quantile(S(1,:), lower_quantile);
% q99 = quantile(S(1,:), 0.99);
S = S(:,(S(1,:) > q1));

% tollerance = 1e-3;
N = size(S,1);

ll = 0;
prev_ll = -1;

s1 = S(1,:);
flags = s1~=0;
s1 = s1(flags);


M1 = size(s1,2);
 
S1_sorted = sort(s1, 'descend');
% S_sorted = S;

prior_thres = S1_sorted(1) - (S1_sorted(1) - S1_sorted(M1)) * c_range;
prior_thres2 = S1_sorted(1) - (S1_sorted(1) - S1_sorted(M1)) * (c_range + (1-c_range)/2);
nc1 = sum(S1_sorted > prior_thres);
nic = sum(S1_sorted > prior_thres2 & S1_sorted <= prior_thres);

alpha = nc1 / M1;
beta = nic / M1;


[u_c, sigma_c] = normfit(S1_sorted(1:nc1));
[u_ic, sigma_ic] = normfit(S1_sorted(nc1:nic));
[u_i1, sigma_i1] = normfit(S1_sorted(nic:M1));

lambda_c  = sl1;
lambda_ic = sl2;
lambda_i1 = sl3;

x_range = S1_sorted(M1):(S1_sorted(1)-S1_sorted(M1))/100:S1_sorted(1);
param_c_ic_check_func = make_param_check_function( ...
    x_range, 1, 1, ...
    u_c, sigma_c, lambda_c, ...
    0, sigma_c, lambda_ic, ...
    constraints);

lambda_ic = param_bin_search(1e-7*sign(lambda_ic), lambda_ic, ...
    @(x, fuzzy) param_c_ic_check_func('lambda_2', x, fuzzy));
sigma_ic = param_bin_search(sigma_c, sigma_ic, ...
    @(x, fuzzy) param_c_ic_check_func('sigma_2', x, fuzzy));
u_ic = param_bin_search(0, u_ic, ...
    @(x, fuzzy) param_c_ic_check_func('u_2', x, fuzzy));

param_ic_i1_check_func = make_param_check_function( ...
    x_range, 1, 1, ...
    u_ic, sigma_ic, lambda_ic, ...
    0, sigma_ic, lambda_ic, ...
    constraints);

lambda_i1 = param_bin_search(lambda_ic, lambda_i1, ...
    @(x, fuzzy) param_ic_i1_check_func('lambda_2', x, fuzzy));
sigma_i1 = param_bin_search(sigma_ic, sigma_i1, ...
    @(x, fuzzy) param_ic_i1_check_func('sigma_2', x, fuzzy));
u_i1 = param_bin_search(0, u_i1, ...
    @(x, fuzzy) param_ic_i1_check_func('u_2', x, fuzzy));

SIM_SIZE = 1000;

w1c  = alpha;
w1ic = beta;
w1i1 = 1 - alpha - beta

ws = [w1c, w1ic, w1i1];

figure('Position', [10,10,2000,500]);
while abs(ll - prev_ll) > tolerance
    prev_ll = ll;

    u_c_new  = u_c;
    u_ic_new = u_ic;
    u_i1_new = u_i1;

    sigma_c_new  = sigma_c;
    sigma_ic_new = sigma_ic;
    sigma_i1_new = sigma_i1;

    lambda_c_new  = lambda_c;
    lambda_ic_new = lambda_ic;
    lambda_i1_new = lambda_i1;

    plot_dist_1i_xl_fn(S, '_s1_c3_xl', ws, ...
        u_c, sigma_c, lambda_c, ...
        u_ic, sigma_ic, lambda_ic, ...
        u_i1, sigma_i1, lambda_i1);
    pause(.001);

    ll = func_ll1_3i_xl(s1, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new)
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
    
    [Vc1 , Wc1 ] = trunc_norm_moments(delta_c  / sigma_c  * (s1-u_c ), sqrt(1-delta_c ^2));
    [Vic1, Wic1] = trunc_norm_moments(delta_ic / sigma_ic * (s1-u_ic), sqrt(1-delta_ic^2));
    [Vi11, Wi11] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));
    
    Rs1 = rs_skewnorm(s1, ...
        [w1c, w1ic, w1i1], ...
        [u_c, u_ic, u_i1], ...
        [sigma_c, sigma_ic, sigma_i1], ...
        [lambda_c, lambda_ic, lambda_i1]);
    
    Rs1c  = Rs1(1,:);
    Rs1ic = Rs1(2,:);
    Rs1i1 = Rs1(3,:);
    
    sum_Rs1c  = sum(Rs1c);
    sum_Rs1ic = sum(Rs1ic);
    sum_Rs1i1 = sum(Rs1i1);
    
    w1c  = sum_Rs1c  / M1;
    w1ic = sum_Rs1ic / M1;
    w1i1 = sum_Rs1i1 / M1;
    ws = [w1c, w1ic, w1i1];
    
    
    u_c_new  = (sum( Rs1c  .* (s1 - Vc1  * Delta_c ) )) / (sum_Rs1c);
    u_ic_new = (sum( Rs1ic .* (s1 - Vic1 * Delta_ic) )) / (sum_Rs1ic);
    u_i1_new = (sum( Rs1i1 .* (s1 - Vi11 * Delta_i1) )) / (sum_Rs1i1);

    u_c_new = param_bin_search(u_c, u_c_new, ...
        @(x, fuzzy) param_c_ic_check_func('u_1', x, fuzzy));
    u_ic_new = param_bin_search(u_ic, u_ic_new, ...
        @(x, fuzzy) param_c_ic_check_func('u_2', x, fuzzy));
    % u_ic_new = param_bin_search(u_ic, u_ic_new, ...
    %     @(x, fuzzy) param_ic_i1_check_func('u_1', x, fuzzy));
    u_i1_new = param_bin_search(u_i1, u_i1_new, ...
        @(x, fuzzy) param_ic_i1_check_func('u_2', x, fuzzy));

    ll = func_ll1_3i_xl(s1, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new)

    Delta_c_new  = (sum( Rs1c  .* Vc1  .* (s1 - u_c_new ) ) ...
                    ) / (sum(Rs1c  .* Wc1 ));
    Delta_ic_new = (sum( Rs1ic .* Vic1 .* (s1 - u_ic_new) ) ...
                    ) / (sum(Rs1ic .* Wic1));
    Delta_i1_new = (sum( Rs1i1 .* Vi11 .* (s1 - u_i1_new) ) ...
                    ) / (sum(Rs1i1 .* Wi11));

    Gamma_c_new  = ( sum( Rs1c .* ((s1 - u_c_new).^2 - ...
                                  2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + ...
                                  Wc1 * Delta_c_new^2) ) ...
                     ) / (sum_Rs1c);
    Gamma_ic_new = ( sum( Rs1ic .* ((s1 - u_ic_new).^2 - ...
                                    2 * Vic1 .* (s1 - u_ic_new) * Delta_ic_new + ...
                                    Wic1 * Delta_ic_new^2) ) ...
                     ) / (sum_Rs1ic);
    Gamma_i1_new = ( sum( Rs1i1 .* ((s1 - u_i1_new).^2 - ...
                                    2 * Vi11 .* (s1 - u_i1_new) * Delta_i1_new + ...
                                    Wi11 * Delta_i1_new^2) ) ...
                     ) / (sum_Rs1i1);

    lambda_c_new = sign(Delta_c_new) * sqrt(Delta_c_new^2 / Gamma_c_new);
    lambda_ic_new = sign(Delta_ic_new) * sqrt(Delta_ic_new^2 / Gamma_ic_new);
    lambda_i1_new = sign(Delta_i1_new) * sqrt(Delta_i1_new^2 / Gamma_i1_new);
    
    sigma_c_new = sqrt(Gamma_c_new + Delta_c_new^2);
    sigma_ic_new = sqrt(Gamma_ic_new + Delta_ic_new^2);
    sigma_i1_new = sqrt(Gamma_i1_new + Delta_i1_new^2); 
    
    ll = func_ll1_3i_xl(s1, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new)

    sigma_c_new = param_bin_search(sigma_c, sigma_c_new, ...
        @(x, fuzzy) param_c_ic_check_func('sigma_1', x, fuzzy));
    sigma_ic_new = param_bin_search(sigma_ic, sigma_ic_new, ...
        @(x, fuzzy) param_c_ic_check_func('sigma_2', x, fuzzy));
    % sigma_ic_new = param_bin_search(sigma_ic, sigma_ic_new, ...
    %     @(x, fuzzy) param_ic_i1_check_func('sigma_1', x, fuzzy));
    sigma_i1_new = param_bin_search(sigma_i1, sigma_i1_new, ...
        @(x, fuzzy) param_ic_i1_check_func('sigma_2', x, fuzzy));

    lambda_c_new = param_bin_search(lambda_c, lambda_c_new, ...
        @(x, fuzzy) param_c_ic_check_func('lambda_1', x, fuzzy));
    lambda_ic_new = param_bin_search(lambda_ic, lambda_ic_new, ...
        @(x, fuzzy) param_c_ic_check_func('lambda_2', x, fuzzy));
    % lambda_ic_new = param_bin_search(lambda_ic, lambda_ic_new, ...
    %     @(x, fuzzy) param_ic_i1_check_func('lambda_1', x, fuzzy));
    lambda_i1_new = param_bin_search(lambda_i1, lambda_i1_new, ...
        @(x, fuzzy) param_ic_i1_check_func('lambda_2', x, fuzzy));

    
    ll = func_ll1_3i_xl(s1, ws, ...
    u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
    u_i1_new, sigma_i1_new, lambda_i1_new)


%     alpha_new = w1 + w2 + w3 + w4;
%     beta_new = w1 + w3 + w5 + w6;
    

    
%     alpha = alpha_new;
%     beta = beta_new;
    u_c = u_c_new;
    u_ic = u_ic_new;
    u_i1 = u_i1_new;

	sigma_c = sigma_c_new;
	sigma_ic = sigma_ic_new;
	sigma_i1 = sigma_i1_new;

	lambda_c = lambda_c_new;
	lambda_ic = lambda_ic_new;
	lambda_i1 = lambda_i1_new;
    
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
