function [alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i, sigma_i, lambda_i] = EM1_3_xl(S,sl1,sl2,sl3)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + beta*fi + (1-alpha-beta)*fic
tollerance = 1e-8;
% tollerance = 1e-3;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

S1_sorted = sort(S(1,:), 'descend');
% S_sorted = S;

prior_thres = 2;
nc = sum(S1_sorted > prior_thres);
% alpha = 0.1;
alpha = nc / M;
% beta = 0.05;
beta = (1-alpha) * 2/3; % assume the number of II is 2 times of IC

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:M*alpha));
[u_i, sigma_i, lambda_i] = sn_para_est(S1_sorted(M*alpha+1:M));
lambda_c = lambda_c * sl1;
lambda_i = lambda_i * sl2;
% lambda_ic = lambda_i * sl3;
% sigma_c = sigma_i2;
u_ic = (u_c + u_i) / 2;
sigma_ic = sigma_i;
lambda_ic = sl3;

s1 = S(1,:);
s1 = s1(s1~=0);
% s2 = S(2,:);
% s2 = s2(s2~=0);

M1 = size(s1,2);
% M2 = size(s2,2);
% u_c = 5;
% sigma_c = 1;
% lambda_c = 4;
% u_i = 0;
% sigma_i = 1;
% lambda_i = 2;
% alpha = 0.5;

while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

    ll = func_ll1_3(s1, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_ic, sigma_ic, lambda_ic);
    disp(ll);
    disp(ll - prev_ll);
    disp([alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_ic, sigma_ic, lambda_ic]);
    
    plot_dist_xl_fn(S, '_1s3c_xl', alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_ic, sigma_ic, lambda_ic);
    pause(.1);
    
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);
    delta_i = lambda_i / sqrt(1+lambda_i^2);
    Delta_i = sigma_i * delta_i;
    Gamma_i = sigma_i^2 * (1-delta_i^2);
    delta_ic = lambda_ic / sqrt(1+lambda_ic^2);
    Delta_ic = sigma_ic * delta_ic;
    Gamma_ic = sigma_ic^2 * (1-delta_ic^2);

    [Vc, Wc] = trunc_norm_moments(delta_c / sigma_c * (s1-u_c), sqrt(1-delta_c^2));
    [Vic, Wic] = trunc_norm_moments(delta_ic / sigma_ic * (s1-u_ic), sqrt(1-delta_ic^2));
    [Vi, Wi] = trunc_norm_moments(delta_i / sigma_i * (s1-u_i), sqrt(1-delta_i^2));

    Rs1 = rs1_3(s1, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_ic, sigma_ic, lambda_ic);
    Rsc = Rs1(1,:);
    Rsi = Rs1(2,:);
    Rsic = Rs1(3,:);
    sum_Rsc = sum(Rsc);
    sum_Rsi = sum(Rsi);
    sum_Rsic = sum(Rsic);
    alpha_new = 1 / M1 * sum_Rsc;
    beta_new = 1 / M1 * sum_Rsi;

    u_c_new = (sum( Rsc .* (s1 - Vc * Delta_c))) / (sum_Rsc);
    u_i_new = (sum( Rsi .* (s1 - Vi * Delta_i) )) / (sum_Rsi);
    u_ic_new = (sum( Rsic .* (s1 - Vic * Delta_ic) )) / (sum_Rsic);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)


    Delta_c_new = (sum(Rsc .* Vc .* (s1 - u_c_new))) / (sum(Rsc .* Wc));
    Delta_i_new = (sum(Rsi .* Vi .* (s1 - u_i_new))) / (sum(Rsi .* Wi));
    Delta_ic_new = (sum(Rsic .* Vic .* (s1 - u_ic_new))) / (sum(Rsic .* Wic));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)

    Gamma_c_new = ( sum( Rsc .* ((s1 - u_c_new).^2 - 2 * Vc .* (s1 - u_c_new) * Delta_c_new + Wc * Delta_c_new^2) ) ) / (sum_Rsc);
    Gamma_i_new = ( sum( Rsi .* ((s1 - u_i_new).^2 - 2 * Vi .* (s1 - u_i_new) * Delta_i_new + Wi * Delta_i_new^2) ) ) / (sum_Rsi);
    Gamma_ic_new = ( sum( Rsic .* ((s1 - u_ic_new).^2 - 2 * Vic .* (s1 - u_ic_new) * Delta_ic_new + Wic * Delta_ic_new^2) ) ) / (sum_Rsic);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)

    alpha = alpha_new;
    beta = beta_new;
    
    u_c = u_c_new;
    Delta_c = Delta_c_new;
    Gamma_c = Gamma_c_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    
    u_i = u_i_new;
    Delta_i = Delta_i_new;
    Gamma_i = Gamma_i_new;
    lambda_i = sign(Delta_i) * sqrt(Delta_i^2 / Gamma_i);
    sigma_i = sqrt(Gamma_i + Delta_i^2);
    
    u_ic = u_ic_new;
    Delta_ic = Delta_ic_new;
    Gamma_ic = Gamma_ic_new;
    lambda_ic = sign(Delta_ic) * sqrt(Delta_ic^2 / Gamma_ic);
    sigma_ic = sqrt(Gamma_ic + Delta_ic^2);
    
end

end
