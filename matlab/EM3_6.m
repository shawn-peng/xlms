function [alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3] = EM3_6(S,sl1,sl2,sl3,sl4)
%EM The EM algorithm to estimate parameters
%   Detailed explanation goes here
% tollerance = 1e-8;
tollerance = 1e-8;
% [N, M] = size(S);
N = size(S,1);

ll = 0;
prev_ll = -1;

run('weight_para_solver.m')

s1 = S(1,:);
s1 = s1(s1~=0);
s2 = S(2,:);
s2 = s2(s2~=0);
s3 = S(3,:);
s3 = s3(s3~=0);

M1 = size(s1,2);
M2 = size(s2,2);
M3 = size(s3,2);

% S1_sorted = sort(S(1,:), 'descend');
% S2_sorted = sort(S(2,:), 'descend');
% S3_sorted = sort(S(3,:), 'descend');
S1_sorted = sort(s1, 'descend');
S2_sorted = sort(s2, 'descend');
S3_sorted = sort(s3, 'descend');
% S_sorted = S;

prior_thres = 20;
nc1 = sum(S1_sorted > prior_thres);
nc2 = sum(S2_sorted > prior_thres);
% alpha = 0.1;
alpha = nc1 / M1;
% beta = 0.05;
beta = nc2 / M2;

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:int32(M1*alpha)));
[u_i1, sigma_i1, lambda_i1] = sn_para_est(S1_sorted(int32(M1*alpha):M1));
[u_i2, sigma_i2, lambda_i2] = sn_para_est(S2_sorted(int32(M2*alpha):M2));
[u_i3, sigma_i3, lambda_i3] = sn_para_est(S3_sorted(int32(M3*alpha):M3));
lambda_c = abs(lambda_c) * sl1;
lambda_i1 = abs(lambda_i1) * sl2;
lambda_i2 = abs(lambda_i2) * sl3;
lambda_i3 = abs(lambda_i3) * sl4;
% lambda_c = sl1;
% lambda_i1 = sl2;
% lambda_i2 = sl3;
% lambda_i3 = sl4;
% sigma_c = sigma_i2;
% u_i1 = u_i2;
% sigma_i1 = sigma_i2;
% lambda_i1 = lambda_i2;

sym_M1 = M1;
sym_M2 = M2;
sym_M3 = M3;
% u_c = 5;
% sigma_c = 1;
% lambda_c = 4;
% u_i = 0;
% sigma_i = 1;
% lambda_i = 2;
% alpha = 0.5;

figure('Position', [10,10,2000,1800]);

while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

%     P_c = skew_norm_pdf(S, u_c, sigma_c, lambda_c);
%     P_i = skew_norm_pdf(S, u_i, sigma_i, lambda_i);
    
    % inner sum sum across rows (k = 1, N) producing 1 x M vector
%     ll = sum( sum(log(P_i), 1) + log(alpha * sum(P_c ./ P_i, 1) + (1-alpha) * N) );
    ll = func_ll3_5(s1, s2, s3, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    disp(ll);
    disp(ll - prev_ll);
    disp([alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3]);
    
    plot_dist_fn(S, '_3s4ci', alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    pause(.001);
    
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);
    delta_i2 = lambda_i2 / sqrt(1+lambda_i2^2);
    Delta_i2 = sigma_i2 * delta_i2;
    Gamma_i2 = sigma_i2^2 * (1-delta_i2^2);
    delta_i3 = lambda_i3 / sqrt(1+lambda_i3^2);
    Delta_i3 = sigma_i3 * delta_i3;
    Gamma_i3 = sigma_i3^2 * (1-delta_i3^2);


    [Vc1, Wc1] = trunc_norm_moments(delta_c / sigma_c * (s1-u_c), sqrt(1-delta_c^2));
    [Vc2, Wc2] = trunc_norm_moments(delta_c / sigma_c * (s2-u_c), sqrt(1-delta_c^2));
    [Vi11, Wi11] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));
    [Vi21, Wi21] = trunc_norm_moments(delta_i1 / sigma_i1 * (s2-u_i1), sqrt(1-delta_i1^2));
    [Vi31, Wi31] = trunc_norm_moments(delta_i1 / sigma_i1 * (s3-u_i1), sqrt(1-delta_i1^2));
    [Vi22, Wi22] = trunc_norm_moments(delta_i2 / sigma_i2 * (s2-u_i2), sqrt(1-delta_i2^2));
    [Vi32, Wi32] = trunc_norm_moments(delta_i2 / sigma_i2 * (s3-u_i2), sqrt(1-delta_i2^2));
    [Vi33, Wi33] = trunc_norm_moments(delta_i3 / sigma_i3 * (s3-u_i3), sqrt(1-delta_i3^2));
   
%     S1 = S(1,:);
%     S2 = S(2,:);

%     Rs1c = Rsn1(1,:);
%     Rs1i1 = Rsn1(2,:);
%     Rs2c = Rsn2(1,:);
%     Rs2i1 = Rsn2(2,:);
%     Rs2i2 = Rsn2(3,:);
%     Rs3i2 = Rsn3(1,:);
%     Rs3i3 = Rsn3(2,:);
    [Rs1c, Rs1i1, Rs2c, Rs2i1, Rs2i2, Rs3i2, Rs3i3] = ...
        rsdn(S, [alpha, beta, 1-alpha-beta], ...
            [u_c, u_i1, u_i2, u_i3], ...
            [sigma_c, sigma_i1, sigma_i2, sigma_i3], ...
            [lambda_c, lambda_i1, lambda_i2, lambda_i3]);
    
    sum_Rs1c = sum(Rs1c);
    sum_Rs1i1 = sum(Rs1i1);
    sum_Rs2c = sum(Rs2c);
    sum_Rs2i1 = sum(Rs2i1);
    sum_Rs2i2 = sum(Rs2i2);
    sum_Rs3i2 = sum(Rs3i2);
    sum_Rs3i3 = sum(Rs3i3);
    
%     assert(sum_Rs2c + sum_Rs2i1 + sum_Rs2i2 - M2 < 1e-4)

%     alpha_new = (sum(Rs1c) + sum(Rs2i1)) / (M1 + M2);
%     beta_new = sum(Rs2c) * (1 + (M1 - sum_Rs1c) / (M2 - sum_Rs2i1)) / (M1 + M2);
    C1 = sum_Rs1c;
    I11 = sum_Rs1i1;
    C2 = sum_Rs2c;
    I21 = sum_Rs2i1;
    I22 = sum_Rs2i2;
    I32 = sum_Rs3i2;
    I33 = sum_Rs3i3;
%     alpha_new = double(subs(weights_eqs.sym_alpha));
%     beta_new = double(subs(weights_eqs.sym_beta));
%     flags = alpha_new > 0 & beta_new > 0;
%     alpha_new = alpha_new(flags);
%     beta_new = beta_new(flags);
    alpha_new = C1 / M1;
    beta_new = C2 / M2;

    u_c_new = (sum( Rs1c .* (s1 - Vc1 * Delta_c)) + sum(Rs2c .* (s2 - Vc2 * Delta_c))) / (sum_Rs1c + sum_Rs2c);
    u_i1_new = (sum( Rs1i1 .* (s1 - Vi11 * Delta_i1) ) + sum( Rs2i1 .* (s2 - Vi21 * Delta_i1))) / (sum_Rs1i1 + sum_Rs2i1);
    u_i2_new = (sum( Rs2i2 .* (s2 - Vi22 * Delta_i2) ) + sum( Rs3i2 .* (s3 - Vi32 * Delta_i2) )) / (sum_Rs2i2 + sum_Rs3i2);
    u_i3_new = (sum( Rs3i3 .* (s3 - Vi33 * Delta_i3) )) / (sum_Rs3i3);
%     ll = func_ll2(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)


    Delta_c_new = (sum(Rs1c .* Vc1 .* (s1 - u_c_new)) + sum(Rs2c .* Vc2 .* (s2 - u_c_new))) / (sum(Rs1c .* Wc1) + sum(Rs2c .* Wc2));
    Delta_i1_new = (sum(Rs1i1 .* Vi11 .* (s1 - u_i1_new)) + sum(Rs2i1 .* Vi21 .* (s2 - u_i1_new))) / (sum(Rs1i1 .* Wi11) + sum(Rs2i1 .* Wi21));
    Delta_i2_new = (sum(Rs2i2 .* Vi22 .* (s2 - u_i2_new)) + sum(Rs3i2 .* Vi32 .* (s3 - u_i2_new)) ) / (sum(Rs2i2 .* Wi22) + sum(Rs3i2 .* Wi32));
    Delta_i3_new = (sum(Rs3i3 .* Vi33 .* (s3 - u_i3_new))) / (sum(Rs3i3 .* Wi33));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)
    
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i_new, Gamma_i)

    Gamma_c_new = ( sum( Rs1c .* ((s1 - u_c_new).^2 - 2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + Wc1 * Delta_c_new^2) ) ...
        + sum( Rs2c .* ((s2 - u_c_new).^2 - 2 * Vc2 .* (s2 - u_c_new) * Delta_c_new + Wc2 * Delta_c_new^2) ) ) ...
        / (sum_Rs1c + sum_Rs2c);
    Gamma_i1_new = ( sum( Rs1i1 .* ((s1 - u_i1_new).^2 - 2 * Vi11 .* (s1 - u_i1_new) * Delta_i1_new + Wi11 * Delta_i1_new^2) ) ...
        + sum( Rs2i1 .* ((s2 - u_i1_new).^2 - 2 * Vi21 .* (s2 - u_i1_new) * Delta_i1_new + Wi21 * Delta_i1_new^2) ) ) ...
        / (sum_Rs1i1 + sum_Rs2i1);
    Gamma_i2_new = ( sum( Rs2i2 .* ((s2 - u_i2_new).^2 - 2 * Vi22 .* (s2 - u_i2_new) * Delta_i2_new + Wi22 * Delta_i2_new^2) ) ...
        + sum( Rs3i2 .* ((s3 - u_i2_new).^2 - 2 * Vi32 .* (s3 - u_i2_new) * Delta_i2_new + Wi32 * Delta_i2_new^2) ) ) ...
        / (sum_Rs2i2 + sum_Rs3i2);
    Gamma_i3_new = ( sum( Rs3i3 .* ((s3 - u_i3_new).^2 - 2 * Vi33 .* (s3 - u_i3_new) * Delta_i3_new + Wi33 * Delta_i3_new^2) ) ) ...
        / sum_Rs3i3;
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)

%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i_new)

    alpha = alpha_new;
    beta = beta_new;
    
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
    
    u_i2 = u_i2_new;
    Delta_i2 = Delta_i2_new;
    Gamma_i2 = Gamma_i2_new;
    lambda_i2 = sign(Delta_i2) * sqrt(Delta_i2^2 / Gamma_i2);
    sigma_i2 = sqrt(Gamma_i2 + Delta_i2^2);
    
    u_i3 = u_i3_new;
    Delta_i3 = Delta_i3_new;
    Gamma_i3 = Gamma_i3_new;
    lambda_i3 = sign(Delta_i3) * sqrt(Delta_i3^2 / Gamma_i3);
    sigma_i3 = sqrt(Gamma_i3 + Delta_i3^2);
    
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
