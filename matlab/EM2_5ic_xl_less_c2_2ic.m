function [params, ll, ll1] = EM2_5ic_xl_less_c2_2ic(S,sl1,sl2,sl3,sl4,sl5,varargin)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
%   f2 = alpha*fi1 + beta*fc + (1-alpha-beta)*fi2

parser = inputParser;
addParameter(parser, 'tolerance', 1e-7);
addParameter(parser, 'lower_quantile', 0.01);
addParameter(parser, 'c_range', 1/2);

parse(parser, varargin{:});
tolerance      = parser.Results.tolerance;
lower_quantile  = parser.Results.lower_quantile;
c_range         = parser.Results.c_range;


S = S(:, S(2,:)~=0);
q1 = quantile(S(1,:), lower_quantile);
S = S(:, (S(1,:) > q1));

N = size(S,1);

ll = 0;
prev_ll = -1;

s1 = S(1,:);
s1 = s1(s1~=0);
s2 = S(2,:);
s2 = s2(s2~=0);

sorted_s2 = sort(s2);

M1 = size(s1,2);
M2 = size(s2,2);


[alpha, beta, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1] = ...
    EM1_3ic_xl(S, sl1, sl2, sl3, ...
        'tolerance', tolerance, ...
        'lower_quantile', lower_quantile, ...
        'c_range', c_range);
close;

w1c = alpha;
w1ic = beta;
w1i1 = 1 - alpha - beta;


Rs1 = rsn(s1, ...
        [w1c, w1ic, w1i1], ...
        [u_c, u_ic, u_i1], ...
        [sigma_c, sigma_ic, sigma_i1], ...
        [lambda_c, lambda_ic, lambda_i1]);

Rs1c = Rs1(1, :);
Rs1ic = Rs1(2, :);
Rs1i1 = Rs1(3, :);

% a rough estimation
w2c = 0.01;
w2ic = (w1i1+w1ic)/2;
w2ic2 = w1c / 2;
w2i1 = (w1c+w1ic)/2;
w2i2 = 1 - w2c - w2ic - w2ic2 - w2i1;

ws = [w1c, w1ic, w1i1, w2c, w2ic, w2ic2, w2i1, w2i2];

assert(w2c+w2ic+w2ic2+w2i1+w2i2 - 1.0 <= 1e-4);
assert(all(ws >= 0));

lambda_i2 = sl4;
lambda_ic2 = sl5;
[u_i2, sigma_i2, lambda_i2] = fit_skewnorm_weighted(s2, Rs1i1, u_i1, sigma_i1, lambda_i2);
[u_ic2, sigma_ic2, lambda_ic2] = fit_skewnorm_weighted(s2, Rs1c, u_c, sigma_c, lambda_ic2);


x_range = sorted_s2(M1) : (sorted_s2(1)-sorted_s2(M1))/100 : sorted_s2(1);
param_c_ic_check_func = make_param_check_function_s2combined( ...
    x_range, 1, 1, w2c, w2ic, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic);
param_ic_i1_check_func = make_param_check_function( ...
    x_range, 1, 1, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1);
param_i1_i2_check_func = make_param_check_function( ...
    x_range, 1, 1, ...
    u_i1, sigma_i1, lambda_i1, ...
    u_i2, sigma_i2, lambda_i2);

if ~param_c_ic_check_func('sigma_1', sigma_c, 1e-7)
    new_w2c = param_bin_search(1e-12, w2c, ...
        @(x, fuzzy) param_c_ic_check_func('w_21', x, fuzzy));
    rw = new_w2c / (new_w2c + w2ic);
    restw = w2c - new_w2c;
    w2c = new_w2c + rw * restw;
    w2ic = w2ic + (1-rw) * restw;
    param_c_ic_check_func('w_21', w2c, 1e-7);
    assert(param_c_ic_check_func('w_22', w2ic, 1e-7));
end

if ~param_ic_i1_check_func('sigma_2', sigma_i1, 1e-7)
    
    param_ic_i1_check_func = make_param_check_function( ...
        x_range, 1, 1, ...
        u_ic, sigma_ic, lambda_ic, ...
        u_i1, sigma_ic, lambda_i1);
    disp('Searching lambda_i1');
    
    lambda_i1 = param_bin_search( ...
        1e-7*lambda_i1, lambda_i1, ...
        @(x, fuzzy) param_ic_i1_check_func('lambda_2', x, fuzzy));
    sigma_i1 = param_bin_search( ...
        sigma_ic, sigma_i1, ...
        @(x, fuzzy) param_ic_i1_check_func('sigma_2', x, fuzzy));

end

ws = [w1c, w1ic, w1i1, w2c, w2ic, w2ic2, w2i1, w2i2];
assert(w2c+w2ic+w2ic2+w2i1+w2i2 - 1.0 <= 1e-4);


if ~param_i1_i2_check_func('sigma_2', sigma_i2, 1e-7)

    param_i1_i2_check_func = make_param_check_function( ...
        x_range, 1, 1, ...
        u_i1, sigma_i1, 1e-7*sign(lambda_i1), ...
        -10, sigma_i1, 1e-7*sign(lambda_i2));
    

    u_i2 = param_bin_search( ...
        -10, u_i2, ...
        @(x, fuzzy) param_i1_i2_check_func('u_2', x, fuzzy));
    if ~param_i1_i2_check_func('sigma_2', sigma_i1, 1e-10)
        return
    end

    disp('Searching sigma_i2');
    sigma_i2 = param_bin_search( ...
        sigma_i1, sigma_i2, ...
        @(x, fuzzy) param_i1_i2_check_func('sigma_2', x, fuzzy));
    disp('Searching lambda_i2');
    lambda_i2 = param_bin_search( ...
        1e-7*lambda_i2, lambda_i2, ...
        @(x, fuzzy) param_i1_i2_check_func('lambda_2', x, fuzzy));
end


constrain_violates = 0;

% variables
syms sym_w1c sym_w1ic sym_w1i1;
syms sym_w2c sym_w2ic sym_w2ic2 sym_w2i1 sym_w2i2;
syms sym_eta1 sym_eta2;
syms sym_eta2c sym_eta2ic sym_eta2ic2 sym_eta2i1 sym_eta2i2;
% values
syms sym_R1c sym_R1ic sym_R1i1;
syms sym_R2c sym_R2ic sym_R2ic2 sym_R2i1 sym_R2i2;

sys_eqs = [
    sym_R1c   / sym_w1c   + sym_eta1 + sym_eta2ic + sym_eta2ic2 + sym_eta2i1 == 0;
    sym_R1ic  / sym_w1ic  + sym_eta1 + sym_eta2c  + sym_eta2ic  + sym_eta2i1 == 0;
    sym_R1i1  / sym_w1i1  + sym_eta1 + sym_eta2c  + sym_eta2ic  + sym_eta2i2 == 0;
    sym_R2c   / sym_w2c   + sym_eta2 - sym_eta2c   == 0;
    sym_R2ic  / sym_w2ic  + sym_eta2 - sym_eta2ic  == 0;
    sym_R2ic2 / sym_w2ic2 + sym_eta2 - sym_eta2ic2 == 0;
    sym_R2i1  / sym_w2i1  + sym_eta2 - sym_eta2i1  == 0;
    sym_R2i2  / sym_w2i2  + sym_eta2 - sym_eta2i2  == 0;
    sym_w1c   + sym_w1ic  + sym_w1i1 == 1;
    sym_w2c   + sym_w2ic  + sym_w2ic2 + sym_w2i1 + sym_w2i2 == 1;
];

figure('Position', [10,10,2000,500]);
while abs(ll - prev_ll) > tolerance
    prev_ll = ll;

    u_c_new   = u_c;
    u_ic_new  = u_ic;
    u_ic2_new = u_ic2;
    u_i1_new  = u_i1;
    u_i2_new  = u_i2;

    sigma_c_new   = sigma_c;
    sigma_ic_new  = sigma_ic;
    sigma_ic2_new = sigma_ic2;
    sigma_i1_new  = sigma_i1;
    sigma_i2_new  = sigma_i2;

    lambda_c_new   = lambda_c;
    lambda_ic_new  = lambda_ic;
    lambda_ic2_new = lambda_ic2;
    lambda_i1_new  = lambda_i1;
    lambda_i2_new  = lambda_i2;
    plot_dist_i_xl_ic2_fn(S, '_s2_c5_xl', ws, ...
        u_c,   sigma_c,   lambda_c, ...
        u_ic,  sigma_ic,  lambda_ic, ...
        u_ic2, sigma_ic2, lambda_ic2, ...
        u_i1,  sigma_i1,  lambda_i1, ...
        u_i2,  sigma_i2,  lambda_i2);
    pause(.001);

    ll = func_ll2_5i_xl(s1, s2, ws, ...
        u_c_new,   sigma_c_new,   lambda_c_new, ...
        u_ic_new,  sigma_ic_new,  lambda_ic_new, ...
        u_ic2_new, sigma_ic2_new, lambda_ic2_new, ...
        u_i1_new,  sigma_i1_new,  lambda_i1_new, ...
        u_i2_new,  sigma_i2_new,  lambda_i2_new)
    
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);
    delta_ic = lambda_ic / sqrt(1+lambda_ic^2);
    Delta_ic = sigma_ic * delta_ic;
    Gamma_ic = sigma_ic^2 * (1-delta_ic^2);
    delta_ic2 = lambda_ic2 / sqrt(1+lambda_ic2^2);
    Delta_ic2 = sigma_ic2 * delta_ic2;
    Gamma_ic2 = sigma_ic2^2 * (1-delta_ic2^2);
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);
    delta_i2 = lambda_i2 / sqrt(1+lambda_i2^2);
    Delta_i2 = sigma_i2 * delta_i2;
    Gamma_i2 = sigma_i2^2 * (1-delta_i2^2);

    [Vc1  , Wc1  ] = trunc_norm_moments(delta_c   / sigma_c   * (s1-u_c  ), sqrt(1-delta_c  ^2));
    [Vc2  , Wc2  ] = trunc_norm_moments(delta_c   / sigma_c   * (s2-u_c  ), sqrt(1-delta_c  ^2));
    [Vic1 , Wic1 ] = trunc_norm_moments(delta_ic  / sigma_ic  * (s1-u_ic ), sqrt(1-delta_ic ^2));
    [Vic2 , Wic2 ] = trunc_norm_moments(delta_ic  / sigma_ic  * (s2-u_ic ), sqrt(1-delta_ic ^2));
    [Vic22, Wic22] = trunc_norm_moments(delta_ic2 / sigma_ic2 * (s2-u_ic2), sqrt(1-delta_ic2^2));
    [Vi11 , Wi11 ] = trunc_norm_moments(delta_i1  / sigma_i1  * (s1-u_i1 ), sqrt(1-delta_i1 ^2));
    [Vi12 , Wi12 ] = trunc_norm_moments(delta_i1  / sigma_i1  * (s2-u_i1 ), sqrt(1-delta_i1 ^2));
    [Vi22 , Wi22 ] = trunc_norm_moments(delta_i2  / sigma_i2  * (s2-u_i2 ), sqrt(1-delta_i2 ^2));

    Rs1 = rsn(s1, ...
        [w1c, w1ic, w1i1], ...
        [u_c, u_ic, u_i1], [sigma_c, sigma_ic, sigma_i1], ...
        [lambda_c, lambda_ic, lambda_i1]);
    Rs2 = rsn(s2, ...
        [w2c, w2ic, w2ic2, w2i1, w2i2], ...
        [u_c, u_ic, u_ic2, u_i1, u_i2], ...
        [sigma_c, sigma_ic, sigma_ic2, sigma_i1, sigma_i2], ...
        [lambda_c, lambda_ic, lambda_ic2, lambda_i1, lambda_i2]);
    
    Rs1c  = Rs1(1,:);
    Rs1ic = Rs1(2,:);
    Rs1i1 = Rs1(3,:);
    
    Rs2c   = Rs2(1,:);
    Rs2ic  = Rs2(2,:);
    Rs2ic2 = Rs2(3,:);
    Rs2i1  = Rs2(4,:);
    Rs2i2  = Rs2(5,:);
    
    sum_Rs1c  = sum(Rs1c);
    sum_Rs1ic = sum(Rs1ic);
    sum_Rs1i1 = sum(Rs1i1);
    
    sum_Rs2c   = sum(Rs2c);
    sum_Rs2ic  = sum(Rs2ic);
    sum_Rs2ic2 = sum(Rs2ic2);
    sum_Rs2i1  = sum(Rs2i1);
    sum_Rs2i2  = sum(Rs2i2);
    
    
    w1c   = sum_Rs1c   / M1;
    w1ic  = sum_Rs1ic  / M1;
    w1i1  = sum_Rs1i1  / M1;
    w2c   = sum_Rs2c   / M2;
    w2ic  = sum_Rs2ic  / M2;
    w2ic2 = sum_Rs2ic2 / M2;
    w2i1  = sum_Rs2i1  / M2;
    w2i2  = sum_Rs2i2  / M2;
    
    sub_vars = [sym_R1c, sym_R1ic, sym_R1i1, sym_R2c, sym_R2ic, sym_R2ic2, sym_R2i1, sym_R2i2];
    sub_vals = [sum_Rs1c, sum_Rs1ic, sum_Rs1i1, sum_Rs2c, sum_Rs2ic, sum_Rs2ic2, sum_Rs2i1, sum_Rs2i2];
    
    extra_eqs = sym.empty(0);
    flag = 0;
    
    if w2c <= w1ic + w1i1
        extra_eqs(end+1, 1) = sym_eta2c == 0;
    else
        extra_eqs(end+1, 1) = sym_w2c == sym_w1ic + sym_w1i1;
        flag = 1;
    end
    
    if w2ic <= w1c + w1ic + w1i1
        extra_eqs(end+1, 1) = sym_eta2ic == 0;
    else
        extra_eqs(end+1, 1) = sym_w2ic == sym_w1c + sym_w1ic + sym_w1i1;
        flag = 1;
    end
    
    if w2ic2 <= w1c
        extra_eqs(end+1, 1) = sym_eta2ic2 == 0;
    else
        extra_eqs(end+1, 1) = sym_w2ic2 == sym_w1c;
        flag = 1;
    end
    
    if w2i1 <= w1c + w1ic
        extra_eqs(end+1, 1) = sym_eta2i1 == 0;
    else
        extra_eqs(end+1, 1) = sym_w2i1 == sym_w1c + sym_w1ic;
        flag = 1;
    end
    
    if w2i2 <= w1i1
        extra_eqs(end+1, 1) = sym_eta2i2 == 0;
    else
        extra_eqs(end+1, 1) = sym_w2i2 == sym_w1i1;
        flag = 1;
    end
    
    if flag
        new_sys_eqs = subs([sys_eqs; extra_eqs], sub_vars, sub_vals);
        solution = solve(new_sys_eqs);

        w1c = double(solution.sym_w1c);
        w1ic = double(solution.sym_w1ic);
        w1i1 = double(solution.sym_w1i1);
        w2c = double(solution.sym_w2c);
        w2ic = double(solution.sym_w2ic);
        w2ic2 = double(solution.sym_w2ic2);
        w2i1 = double(solution.sym_w2i1);
        w2i2 = double(solution.sym_w2i2);
    end

    % assert(w2c <= w1ic + w1i1);
    % assert(w2ic <= w1c + w1ic + w1i1);
    % assert(w2ic2 <= w1c);
    % assert(w2i1 <= w1c + w1ic);
    % assert(w2i2 <= w1i1);

    old_ws = ws;

    old_w2c = old_ws(4);
    new_w2c = param_bin_search( ...
        old_w2c, w2c, ...
        @(x, fuzzy) param_c_ic_check_func('w_21', x, fuzzy));
    rw = new_w2c / (new_w2c + w2ic);
    old_sum = w2c + w2ic;
    restw = w2c - new_w2c;
    w2c = new_w2c + rw * restw;
    w2ic = w2ic + (1-rw) * restw;
    new_sum = w2c + w2ic;
    w2ic = w2ic + (old_sum - new_sum);
    param_c_ic_check_func('w_21', w2c, 1e-7);
    assert(param_c_ic_check_func('w_22', w2ic, 1e-7));
    
    ws = [w1c, w1ic, w1i1, w2c, w2ic, w2ic2, w2i1, w2i2];
    
    assert(all(ws >= 0));
    assert(abs(sum(ws) - 2) < 1e-10);
    
    ll = func_ll2_4i_xl(s1, s2, ws, ...
    u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
    u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new)


    u_c_new  = (sum( Rs1c  .* (s1 - Vc1  * Delta_c ) ) + ...
                sum( Rs2c  .* (s2 - Vc2  * Delta_c ) )) / (sum_Rs1c + sum_Rs2c);
    u_ic_new = (sum( Rs1ic .* (s1 - Vic1 * Delta_ic) ) + ...
                sum( Rs2ic .* (s2 - Vic2 * Delta_ic) )) / (sum_Rs1ic + sum_Rs2ic);
    u_i1_new = (sum( Rs1i1 .* (s1 - Vi11 * Delta_i1) ) + ...
                sum( Rs2i1 .* (s2 - Vi12 * Delta_i1) )) / (sum_Rs1i1 + sum_Rs2i1);
    u_i2_new = (sum( Rs2i2 .* (s2 - Vi22 * Delta_i2) )) / (sum_Rs2i2);
    u_ic2_new = (sum( Rs2ic2 .* (s2 - Vic22 * Delta_ic2) )) / (sum_Rs2ic2);
    
    if u_i2_new > u_i1_new
        u_i2_new = u_i1_new;
    end
    

    u_c_new = param_bin_search(u_c, u_c_new, @(x, fuzzy) param_c_ic_check_func('u_1', x, fuzzy));
    u_ic_new = param_bin_search(u_ic, u_ic_new, @(x, fuzzy) param_c_ic_check_func('u_2',  x, fuzzy));
    u_ic_new = param_bin_search(u_ic, u_ic_new, @(x, fuzzy) param_ic_i1_check_func('u_1', x, fuzzy));
    u_i1_new = param_bin_search(u_i1, u_i1_new, @(x, fuzzy) param_ic_i1_check_func('u_2', x, fuzzy));
    u_i1_new = param_bin_search(u_i1, u_i1_new, @(x, fuzzy) param_i1_i2_check_func('u_1', x, fuzzy));
    u_i2_new = param_bin_search(u_i2, u_i2_new, @(x, fuzzy) param_i1_i2_check_func('u_2', x, fuzzy));

    ll = func_ll2_4i_xl(s1, s2, ws, ...
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
    Delta_ic2_new = (sum( Rs2ic2 .* Vic22 .* (s2 - u_ic2_new) ) ...
                    ) / (sum(Rs2ic2 .* Wic22));

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
    Gamma_ic2_new = ( sum( Rs2ic2 .* ((s2 - u_ic2_new).^2 - ...
                                    2 * Vic22 .* (s2 - u_ic2_new) * Delta_ic2_new + ...
                                    Wic22 * Delta_ic2_new^2) ) ...
                     ) / (sum_Rs2ic2);

    lambda_c_new = sign(Delta_c_new) * sqrt(Delta_c_new^2 / Gamma_c_new);
    lambda_ic_new = sign(Delta_ic_new) * sqrt(Delta_ic_new^2 / Gamma_ic_new);
    lambda_i1_new = sign(Delta_i1_new) * sqrt(Delta_i1_new^2 / Gamma_i1_new);
    lambda_i2_new = sign(Delta_i2_new) * sqrt(Delta_i2_new^2 / Gamma_i2_new);
    lambda_ic2_new = sign(Delta_ic2_new) * sqrt(Delta_ic2_new^2 / Gamma_ic2_new);
    
    sigma_c_new = sqrt(Gamma_c_new + Delta_c_new^2);
    sigma_ic_new = sqrt(Gamma_ic_new + Delta_ic_new^2);
    sigma_i1_new = sqrt(Gamma_i1_new + Delta_i1_new^2);
    sigma_i2_new = sqrt(Gamma_i2_new + Delta_i2_new^2);
    sigma_ic2_new = sqrt(Gamma_ic2_new + Delta_ic2_new^2);
    
    ll = func_ll2_4i_xl(s1, s2, ws, ...
        u_c_new, sigma_c_new, lambda_c_new, u_ic_new, sigma_ic_new, lambda_ic_new, ...
        u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new)

    
    sigma_c_new = param_bin_search(sigma_c, sigma_c_new, @(x, fuzzy) param_c_ic_check_func('sigma_1', x, fuzzy));
    sigma_ic_new = param_bin_search(sigma_ic, sigma_ic_new, @(x, fuzzy) param_c_ic_check_func('sigma_2', x, fuzzy));
    sigma_ic_new = param_bin_search(sigma_ic, sigma_ic_new, @(x, fuzzy) param_ic_i1_check_func('sigma_1', x, fuzzy));
    sigma_i1_new = param_bin_search(sigma_i1, sigma_i1_new, @(x, fuzzy) param_ic_i1_check_func('sigma_2', x, fuzzy));
    sigma_i1_new = param_bin_search(sigma_i1, sigma_i1_new, @(x, fuzzy) param_i1_i2_check_func('sigma_1', x, fuzzy));
    sigma_i2_new = param_bin_search(sigma_i2, sigma_i2_new, @(x, fuzzy) param_i1_i2_check_func('sigma_2', x, fuzzy));

    lambda_c_new = param_bin_search(lambda_c, lambda_c_new, @(x, fuzzy) param_c_ic_check_func('lambda_1', x, fuzzy));
    lambda_ic_new = param_bin_search(lambda_ic, lambda_ic_new, @(x, fuzzy) param_c_ic_check_func('lambda_2', x, fuzzy));
    lambda_ic_new = param_bin_search(lambda_ic, lambda_ic_new, @(x, fuzzy) param_ic_i1_check_func('lambda_1', x, fuzzy));
    lambda_i1_new = param_bin_search(lambda_i1, lambda_i1_new, @(x, fuzzy) param_ic_i1_check_func('lambda_2', x, fuzzy));
    lambda_i1_new = param_bin_search(lambda_i1, lambda_i1_new, @(x, fuzzy) param_i1_i2_check_func('lambda_1', x, fuzzy));
    lambda_i2_new = param_bin_search(lambda_i2, lambda_i2_new, @(x, fuzzy) param_i1_i2_check_func('lambda_2', x, fuzzy));
    
    u_c   = u_c_new;
    u_ic  = u_ic_new;
    u_i1  = u_i1_new;
    u_i2  = u_i2_new;
    u_ic2 = u_ic2_new;

	sigma_c = sigma_c_new;
	sigma_ic = sigma_ic_new;
	sigma_i1 = sigma_i1_new;
	sigma_i2 = sigma_i2_new;
	sigma_ic2 = sigma_ic2_new;

	lambda_c = lambda_c_new;
	lambda_ic = lambda_ic_new;
	lambda_i1 = lambda_i1_new;
	lambda_i2 = lambda_i2_new;
	lambda_ic2 = lambda_ic2_new;
    
end

theta_c = pack_skntheta(u_c, sigma_c, lambda_c);
theta_ic = pack_skntheta(u_ic, sigma_ic, lambda_ic);
theta_i1 = pack_skntheta(u_i1, sigma_i1, lambda_i1);
theta_i2 = pack_skntheta(u_i2, sigma_i2, lambda_i2);
theta_ic2 = pack_skntheta(u_ic2, sigma_ic2, lambda_ic2);
params = {ws, [theta_c, theta_ic, theta_i1, theta_i2, theta_ic2]};

ll1 = func_ll_skewnorm_mixture(s1, ...
    [w1c, w1ic, w1i1], ...
    [u_c, u_ic, u_i1], ...
    [sigma_c, sigma_ic, sigma_i1], ...
    [lambda_c, lambda_ic, lambda_i1]);

end

function [u, sigma] = unpack2(theta)
    u = theta(1);
    sigma = theta(2);
end

