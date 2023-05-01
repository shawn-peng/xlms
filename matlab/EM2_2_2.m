function [alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1] = EM2_2_2(S,sl1,sl2)
%EM The EM algorithm to estimate parameters
%   f1 = alpha*fc + (1-alpha)*fi1
%   f2 = alpha*fi1 + (1-alpha)*fi2
% tollerance = 1e-8;
tollerance = 1e-5;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

S1_sorted = sort(S(1,:), 'descend');
S2_sorted = sort(S(2,:), 'descend');
% S_sorted = S;

prior_thres = 30;
nc1 = sum(S1_sorted > prior_thres);
nc2 = sum(S2_sorted > prior_thres);
% alpha = 0.1;
alpha = nc1 / M;

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:M*alpha));
[u_i2, sigma_i2, lambda_i2] = sn_para_est(reshape(S2_sorted, 1, []));
lambda_c = lambda_c * sl1;
lambda_i2 = lambda_i2 * sl2;
% sigma_c = sigma_i2;
u_i1 = u_i2;
sigma_i1 = sigma_i2;
lambda_i1 = lambda_i2;

s1 = S(1,:);
s1 = s1(s1~=0);
s2 = S(2,:);
s2 = s2(s2~=0);

M1 = size(s1,2);
M2 = size(s2,2);
% u_c = 5;
% sigma_c = 1;
% lambda_c = 4;
% u_i = 0;
% sigma_i = 1;
% lambda_i = 2;
% alpha = 0.5;

while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

%     P_c = skew_norm_pdf(S, u_c, sigma_c, lambda_c);
%     P_i = skew_norm_pdf(S, u_i, sigma_i, lambda_i);
    
    % inner sum sum across rows (k = 1, N) producing 1 x M vector
%     ll = sum( sum(log(P_i), 1) + log(alpha * sum(P_c ./ P_i, 1) + (1-alpha) * N) );
    ll = func_ll2_3(s1, s2, alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2);
    disp(ll);
    disp(ll - prev_ll);
    disp([alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1]);
    
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);
    delta_i2 = lambda_i2 / sqrt(1+lambda_i2^2);
    Delta_i2 = sigma_i2 * delta_i2;
    Gamma_i2 = sigma_i2^2 * (1-delta_i2^2);
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);

%     [Vc, Wc] = trunc_norm_moments(delta_c / sigma_c * (S-u_c), sqrt(1-delta_c^2));
%     [Vi, Wi] = trunc_norm_moments(delta_i / sigma_i * (S-u_i), sqrt(1-delta_i^2));
% 
%     Vc1 = Vc(1,:);
%     Vc2 = Vc(2,:);
%     Wc1 = Wc(1,:);
%     Wc2 = Wc(2,:);
%     Vi1 = Vi(1,:);
%     Vi2 = Vi(2,:);
%     Wi1 = Wi(1,:);
%     Wi2 = Wi(2,:);

    [Vc1, Wc1] = trunc_norm_moments(delta_c / sigma_c * (s1-u_c), sqrt(1-delta_c^2));
%     [Vc2, Wc2] = trunc_norm_moments(delta_c / sigma_c * (s2-u_c), sqrt(1-delta_c^2));
    [Vi1, Wi1] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));
    [Vi21, Wi21] = trunc_norm_moments(delta_i1 / sigma_i1 * (s2-u_i1), sqrt(1-delta_i1^2));
    [Vi22, Wi22] = trunc_norm_moments(delta_i2 / sigma_i2 * (s2-u_i2), sqrt(1-delta_i2^2));
   
%     S1 = S(1,:);
%     S2 = S(2,:);

    Rs1 = rs(s1, alpha, u_i1, sigma_i1, lambda_i1, u_c, sigma_c, lambda_c);
    Rs2 = rs(s2, alpha, u_i2, sigma_i2, lambda_i2, u_i1, sigma_i1, lambda_i1); % s2 ~ I1
    sum_Rs1 = sum(Rs1);
    sum_Rs2 = sum(Rs2);
    alpha_new = 1 / M * sum(Rs1);

    u_i1_new = (sum( (1 - Rs1) .* (s1 - Vi1 * Delta_i1) ) + sum( Rs2 .* (s2 - Vi21 * Delta_i1) )) / (M1 - sum_Rs1 + sum_Rs2);
    u_i2_new = (sum( (1 - Rs2) .* (s2 - Vi22 * Delta_i2) )) / (M2 - sum_Rs2);
    u_c_new = (sum( Rs1 .* (s1 - Vc1 * Delta_c))) / (sum_Rs1);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)


    Delta_i1_new = (sum((1 - Rs1) .* Vi1 .* (s1 - u_i1_new)) + sum(Rs2 .* Vi21 .* (s2 - u_i1_new))) / (sum((1 - Rs1) .* Wi1) + sum(Rs2 .* Wi21));
    Delta_i2_new = (sum((1 - Rs2) .* Vi22 .* (s2 - u_i2_new))) / (sum((1 - Rs2) .* Wi22));
    Delta_c_new = (sum(Rs1 .* Vc1 .* (s1 - u_c_new))) / (sum(Rs1 .* Wc1));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)
    
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i_new, Gamma_i)

    Gamma_i1_new = ( sum( (1-Rs1) .* ((s1 - u_i1_new).^2 - 2 * Vi1 .* (s1 - u_i1_new) * Delta_i1_new + Wi1 * Delta_i1_new^2) ) ...
        + sum( (Rs2) .* ((s2 - u_i1_new).^2 - 2 * Vi21 .* (s2 - u_i1_new) * Delta_i1_new + Wi21 * Delta_i1_new^2) )) / (M1 - sum_Rs1 + sum_Rs2);
    Gamma_i2_new = ( sum( (1-Rs2) .* ((s2 - u_i2_new).^2 - 2 * Vi22 .* (s2 - u_i2_new) * Delta_i2_new + Wi22 * Delta_i2_new^2) ) ) ...
        / (M2 - sum_Rs2);
    Gamma_c_new = ( sum( Rs1 .* ((s1 - u_c_new).^2 - 2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + Wc1 * Delta_c_new^2) ) ) / (sum_Rs1);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)

%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i_new)

    alpha = alpha_new;
    u_c = u_c_new;
    u_i1 = u_i1_new;
    u_i2 = u_i2_new;
    Delta_c = Delta_c_new;
    Delta_i1 = Delta_i1_new;
    Delta_i2 = Delta_i2_new;
    Gamma_c = Gamma_c_new;
    Gamma_i1 = Gamma_i1_new;
    Gamma_i2 = Gamma_i2_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_i1 = sign(Delta_i1) * sqrt(Delta_i1^2 / Gamma_i1);
    lambda_i2 = sign(Delta_i2) * sqrt(Delta_i2^2 / Gamma_i2);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
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
