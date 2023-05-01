function [alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3] = EM3_3(S,slc,sl1,sl2,sl3)
%EM The EM algorithm to estimate parameters
%   consider the first and the second score are dependent, they can not be
%   true at the same time
% tollerance = 1e-8;
tollerance = 1e-2;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

S1_sorted = sort(S(1,:), 'descend');
S2_sorted = sort(S(2,:), 'descend');
S3_sorted = sort(S(3,:), 'descend');
% S_sorted = S;

prior_thres = 30;
nc1 = sum(S1_sorted > prior_thres);
nc2 = sum(S2_sorted > prior_thres);
% alpha = 0.1;
% alpha = nc1 / M;
alpha = 0.5;
% beta = 0.05;
beta = (nc2 / M) / (1 - alpha);

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:M*alpha));
[u_i1, sigma_i1, lambda_i1] = sn_para_est(S1_sorted(int32(M*alpha):M));
[u_i2, sigma_i2, lambda_i2] = sn_para_est(reshape(S2_sorted, 1, []));
[u_i3, sigma_i3, lambda_i3] = sn_para_est(reshape(S3_sorted, 1, []));
lambda_c = lambda_c * slc;
lambda_i1 = lambda_i1 * sl1;
lambda_i2 = lambda_i2 * sl2;
% sigma_c = sigma_i2;
% u_i3 = u_i2;
% sigma_i3 = sigma_i2;
% lambda_i3 = abs(lambda_i2) * sl3;

u_i = u_i1;
sigma_i = sigma_i1;
lambda_i = lambda_i1;

s1 = S(1,:);
s1 = s1(s1~=0);
s2 = S(2,:);
s2 = s2(s2~=0);
s3 = S(3,:);
s3 = s3(s3~=0);

M1 = size(s1,2);
M2 = size(s2,2);
M3 = size(s3,2);
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
    ll = func_ll3_3(s1, s2, s3, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
%     disp(ll);
    disp(ll - prev_ll);
    disp([alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3]);
    
    plot_dist_fn(S, '_3s4c', alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
%     pause(.001);
    
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);
    delta_i2 = lambda_i2 / sqrt(1+lambda_i2^2);
    Delta_i2 = sigma_i2 * delta_i2;
    Gamma_i2 = sigma_i2^2 * (1-delta_i2^2);
    delta_i3 = lambda_i3 / sqrt(1+lambda_i3^2);
    Delta_i3 = sigma_i3 * delta_i3;
    Gamma_i3 = sigma_i3^2 * (1-delta_i3^2);
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
    [Vc2, Wc2] = trunc_norm_moments(delta_c / sigma_c * (s2-u_c), sqrt(1-delta_c^2));
    [Vi11, Wi11] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));
    [Vi21, Wi21] = trunc_norm_moments(delta_i1 / sigma_i1 * (s2-u_i1), sqrt(1-delta_i1^2));
    [Vi31, Wi31] = trunc_norm_moments(delta_i1 / sigma_i1 * (s3-u_i1), sqrt(1-delta_i1^2));
    [Vi22, Wi22] = trunc_norm_moments(delta_i2 / sigma_i2 * (s2-u_i2), sqrt(1-delta_i2^2));
    [Vi32, Wi32] = trunc_norm_moments(delta_i2 / sigma_i2 * (s3-u_i2), sqrt(1-delta_i2^2));
    [Vi33, Wi33] = trunc_norm_moments(delta_i3 / sigma_i3 * (s3-u_i3), sqrt(1-delta_i3^2));
   
%     S1 = S(1,:);
%     S2 = S(2,:);

    pc1 = skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
    pi11 = skew_norm_pdf(s1, u_i1, sigma_i1, lambda_i1);
    pc2 = skew_norm_pdf(s2, u_c, sigma_c, lambda_c);
    pi21 = skew_norm_pdf(s2, u_i1, sigma_i1, lambda_i1);
    pi22 = skew_norm_pdf(s2, u_i2, sigma_i2, lambda_i2);
    pi32 = skew_norm_pdf(s3, u_i2, sigma_i2, lambda_i2);
    pi33 = skew_norm_pdf(s3, u_i3, sigma_i3, lambda_i3);
    
    p_total = alpha * pc1 .* pi21 .* pi32 + (1-alpha)*beta * pi11 .* pc2 .* pi32 + (1-alpha)*(1-beta) * pi11 .* pi22 .* pi33;
    
    R1 = alpha * pc1 .* pi21 .* pi32 ./ p_total;
    R2 = (alpha*beta * pc1 .* pi21 .* pi32 + (1-alpha)*beta * pi11 .* pc2 .* pi32) ./ p_total;
    
    sum_R1 = sum(R1);
    sum_R2 = sum(R2);
    
    
    alpha_new = sum_R1 / M;
    beta_new = sum_R2 / M;


    u_c_new = (sum( R1 .* (s1 - Vc1 * Delta_c) ) + sum( R2 .* (s2 - Vc2 * Delta_c) )) / (sum_R1 + sum_R2);
    u_i1_new = (sum( (1 - R1) .* (s1 - Vi11 * Delta_i1) ) ...
        + sum( (R1.*(1 - R2)) .* (s2 - Vi21 * Delta_i1) ) ...
        + sum( (R1.*R2) .* (s3 - Vi31 * Delta_i1) )) / M;
    u_i2_new = (sum( ((1 - R1).*(1 - R2)) .* (s2 - Vi22 * Delta_i2) ) ...
        + sum( (R1.*(1 - R2) + (1 - R1).*R2) .* (s3 - Vi32 * Delta_i2) )) / sum(1 - R1.*R2);
    u_i3_new = (sum( ((1 - R1).*(1 - R2)) .* (s3 - Vi33 * Delta_i3) )) / sum((1 - R1).*(1 - R2));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)


    Delta_c_new = (sum( R1 .* Vc1 .* (s1 - u_c_new) ) ...
        + sum( R2 .* Vc2 .* (s2 - u_c_new) ) ...
        ) / (sum(R1 .* Wc1) + sum(R2 .* Wc2));
    Delta_i1_new = (sum( (1 - R1) .* Vi11 .* (s1 - u_i1_new) ) ...
        + sum( (R1.*(1 - R2)) .* Vi21 .* (s2 - u_i1_new) ) ...
        + sum( (R1.*R2) .* Vi31 .* (s3 - u_i1_new) ) ...
        ) / (sum((R1.*(1 - R2)) .* Wi21) + sum((1 - R1) .* Wi11) + sum((R1.*R2) .* Wi31));
    Delta_i2_new = (sum( (R1 + (1 - R1).*R2) .* Vi32 .* (s3 - u_i2_new) ) ...
        + sum( ((1 - R1).*(1 - R2)) .* Vi22 .* (s2 - u_i2_new) ) ...
        ) / (sum((R1.*(1 - R2) + (1 - R1).*R2) .* Wi32) + sum(((1 - R1).*(1 - R2)) .* Wi22));
    Delta_i3_new = (sum( ((1 - R1).*(1 - R2)) .* Vi33 .* (s3 - u_i3_new) )) / (sum( ((1 - R1).*(1 - R2)) .* Wi33));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)
    
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i_new, Gamma_i)

    Gamma_c_new = ( sum( R1 .* ((s1 - u_c_new).^2 - 2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + Wc1 * Delta_c_new^2) ) ...
        + sum( R2 .* ((s2 - u_c_new).^2 - 2 * Vc2 .* (s2 - u_c_new) * Delta_c_new + Wc2 * Delta_c_new^2) ) ) ...
        / (sum_R1 + sum_R2);
    Gamma_i1_new = ( sum( (1 - R1) .* ((s1 - u_i1_new).^2 - 2 * Vi11 .* (s1 - u_i1_new) * Delta_i1_new + Wi11 * Delta_i1_new^2) ) ...
        + sum( R1 .* ((s2 - u_i1_new).^2 - 2 * Vi21 .* (s2 - u_i1_new) * Delta_i1_new + Wi21 * Delta_i1_new^2) ) ...
        + sum( (R1.*R2) .* ((s3 - u_i1_new).^2 - 2 * Vi31 .* (s3 - u_i1_new) * Delta_i1_new + Wi31 * Delta_i1_new^2) ) ) ...
        / M;
    Gamma_i2_new = ( sum( (1 - R1).*(1 - R2) .* ((s2 - u_i2_new).^2 - 2 * Vi22 .* (s2 - u_i2_new) * Delta_i2_new + Wi22 * Delta_i2_new^2) ) ) ...
        + sum( (R1 + (1 - R1).*R2) .* ((s3 - u_i2_new).^2 - 2 * Vi32 .* (s3 - u_i2_new) * Delta_i2_new + Wi32 * Delta_i2_new^2) ) ...
        / M;
    Gamma_i3_new = ( sum( (1 - R1).*(1 - R2) .* ((s3 - u_i3_new).^2 - 2 * Vi33 .* (s3 - u_i3_new) * Delta_i3_new + Wi33 * Delta_i3_new^2) ) ) ...
        / sum((1 - R1).*(1 - R2));
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)

%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i_new)

    alpha = alpha_new;
    beta = beta_new;
    u_c = u_c_new;
    u_i1 = u_i1_new;
    u_i2 = u_i2_new;
    u_i3 = u_i3_new;
    Delta_c = Delta_c_new;
    Delta_i1 = Delta_i1_new;
    Delta_i2 = Delta_i2_new;
    Delta_i3 = Delta_i3_new;
    Gamma_c = Gamma_c_new;
    Gamma_i1 = Gamma_i1_new;
    Gamma_i2 = Gamma_i2_new;
    Gamma_i3 = Gamma_i3_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_i1 = sign(Delta_i1) * sqrt(Delta_i1^2 / Gamma_i1);
    lambda_i2 = sign(Delta_i2) * sqrt(Delta_i2^2 / Gamma_i2);
    lambda_i3 = sign(Delta_i3) * sqrt(Delta_i3^2 / Gamma_i3);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    sigma_i1 = sqrt(Gamma_i1 + Delta_i1^2);
    sigma_i2 = sqrt(Gamma_i2 + Delta_i2^2);
    sigma_i3 = sqrt(Gamma_i3 + Delta_i3^2);
    
%     break
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
