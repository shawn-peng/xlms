function [alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1] = EM2_1(S,sl1,sl2)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
% tollerance = 1e-8;
tollerance = 1e-8;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

S1_sorted = sort(S(1,:), 'descend');
% S2_sorted = sort(S(2,:), 'descend');
% S_sorted = S;

prior_thres = 20;
nc1 = sum(S1_sorted > prior_thres);
% nc2 = sum(S2_sorted > prior_thres);
% alpha = 0.1;
alpha = nc1 / M;

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:int32(M*alpha)));
[u_i1, sigma_i1, lambda_i1] = sn_para_est(S1_sorted(int32(M*alpha):M));

lambda_c = lambda_c * sl1;
lambda_i1 = lambda_i1 * sl2;

s1 = S(1,:);
s1 = s1(s1~=0);

M1 = size(s1,2);

% figure('Position', [10,10,2000,500]);

while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

    ll = func_ll2_1(s1, alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1);
    disp(ll);
    disp(ll - prev_ll);
    disp([alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1]);
    
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);

    [Vc1, Wc1] = trunc_norm_moments(delta_c / sigma_c * (s1-u_c), sqrt(1-delta_c^2));
    [Vi1, Wi1] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));

    Rs1 = rs(s1, alpha, u_i1, sigma_i1, lambda_i1, u_c, sigma_c, lambda_c);
    sum_Rs1 = sum(Rs1);
    alpha_new = 1 / M * sum(Rs1);

    u_c_new = (sum( Rs1 .* (s1 - Vc1 * Delta_c))) / (sum_Rs1);
    u_i1_new = (sum( (1 - Rs1) .* (s1 - Vi1 * Delta_i1) )) / (M1 - sum_Rs1);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)


    Delta_c_new = (sum(Rs1 .* Vc1 .* (s1 - u_c_new))) / (sum(Rs1 .* Wc1));
    Delta_i1_new = (sum((1 - Rs1) .* Vi1 .* (s1 - u_i1_new))) / (sum((1 - Rs1) .* Wi1));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)

    Gamma_c_new = ( sum( Rs1 .* ((s1 - u_c_new).^2 - 2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + Wc1 * Delta_c_new^2) ) ) / (sum_Rs1);
    Gamma_i1_new = ( sum( (1-Rs1) .* ((s1 - u_i1_new).^2 - 2 * Vi1 .* (s1 - u_i1_new) * Delta_i1_new + Wi1 * Delta_i1_new^2) ) ) / (M1 - sum_Rs1);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)

    alpha = alpha_new;
    
    u_c = u_c_new;
    Delta_c = Delta_c_new;
    Gamma_c = Gamma_c_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    
    u_i1 = u_i1_new;
    Delta_i1 = Delta_i1_new;
    Gamma_i1 = Gamma_i1_new;
    lambda_i1 = sign(Delta_i1) * sqrt(Delta_i1^2 / Gamma_i1);
    sigma_i1 = sqrt(Gamma_i1 + Delta_i1^2);
    
end

end
