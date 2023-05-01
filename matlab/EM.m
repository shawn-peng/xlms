function [alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = EM(S)
%EM The EM algorithm to estimate parameters
%   Detailed explanation goes here
tollerance = 1e-8;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

S_sorted = sort(S, 'descend');
S1_sorted = S_sorted(1,:);
prior_thres = 5;
nc1 = sum(S1_sorted > prior_thres);
% alpha = 0.1;
alpha = nc1 / M;
% beta = 0.05;
% beta = nc2 / M;

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:M*alpha));
% [u_c, sigma_c, lambda_c] = sn_para_est(S_sorted(:,1));
[u_i, sigma_i, lambda_i] = sn_para_est(reshape(S_sorted(:,2:end), 1, []));
sigma_c = sigma_i;
alpha = 0.5;

% u_c = 5;
% sigma_c = 1;
% lambda_c = 4;
% u_i = 0;
% sigma_i = 1;
% lambda_i = 2;
% alpha = 0.5;

while abs(ll - prev_ll) > tollerance
    prev_ll = ll;

    P_c = skew_norm_pdf(S, u_c, sigma_c, lambda_c);
    P_i = skew_norm_pdf(S, u_i, sigma_i, lambda_i);
    
    % inner sum sum across rows (k = 1, N) producing 1 x M vector
    ll = sum( sum(log(P_i), 1) + log(alpha * sum(P_c ./ P_i, 1) + (1-alpha) * N) );
    disp(ll);
    disp(ll - prev_ll);
    
    delta_i = lambda_i / sqrt(1+lambda_i^2);
    Delta_i = sigma_i * delta_i;
    Gamma_i = sigma_i^2 * (1-delta_i^2);
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);

    [Vc, Wc] = trunc_norm_moments(delta_c / sigma_c * (S-u_c), sqrt(1-delta_c^2));
    [Vi, Wi] = trunc_norm_moments(delta_i / sigma_i * (S-u_i), sqrt(1-delta_i^2));


    Rs1 = rs_1(S, alpha, u_i, sigma_i, lambda_i, u_c, sigma_c, lambda_c);
    sum_Rs1 = sum(sum(Rs1));
    alpha_new = 1 / M * sum(sum(Rs1));

    u_c_new = sum(sum(Rs1 .* (S - Vc * Delta_c))) / sum_Rs1;
    u_i_new = sum(sum((1 - Rs1) .* (S - Vi * Delta_i))) / (N*M - sum_Rs1);
%     ll = func_ll(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)


    Delta_c_new = sum(sum(Rs1 .* Vc .* (S - u_c_new))) / sum(sum(Rs1 .* Wc));
%     ll = func_ll(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)
    Delta_i_new = sum(sum((1 - Rs1) .* Vi .* (S - u_i_new))) / sum(sum((1 - Rs1) .* Wi));
    
%     ll = func_ll(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i_new, Gamma_i)

    Gamma_c_new = sum(sum( Rs1 .* ((S - u_c_new).^2 - 2 * Vc .* (S - u_c_new) * Delta_c_new + Wc * Delta_c_new^2) )) / sum_Rs1;
%     ll = func_ll(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)
    Gamma_i_new = sum(sum( (1-Rs1) .* ((S - u_i_new).^2 - 2 * Vi .* (S - u_i_new) * Delta_i_new + Wi * Delta_i_new^2) )) / (N*M - sum_Rs1);

%     ll = func_ll(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i_new)

    alpha = alpha_new;
    u_c = u_c_new;
    u_i = u_i_new;
    Delta_c = Delta_c_new;
    Delta_i = Delta_i_new;
    Gamma_c = Gamma_c_new;
    Gamma_i = Gamma_i_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_i = sign(Delta_i) * sqrt(Delta_i^2 / Gamma_i);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    sigma_i = sqrt(Gamma_i + Delta_i^2);
    
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

function ll = func_ll(S, alpha, u_c, Delta_c, Gamma_c, u_i, Delta_i, Gamma_i)

    [N, M] = size(S);
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_i = sign(Delta_i) * sqrt(Delta_i^2 / Gamma_i);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    sigma_i = sqrt(Gamma_i + Delta_i^2);
    
    P_c = skew_norm_pdf(S, u_c, sigma_c, lambda_c);
    P_i = skew_norm_pdf(S, u_i, sigma_i, lambda_i);
    
    % inner sum sum across rows (k = 1, N) producing 1 x M vector
    ll = sum( sum(log(P_i), 1) + log(alpha * sum(P_c ./ P_i, 1) + (1-alpha) * N) );
end
