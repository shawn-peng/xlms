function [u_new, sigma_new, lambda_new] = fit_skewnorm_weighted_one_step(s, Rs, u, sigma, lambda)
    tolerance = 1e-12;
    sum_diff = inf;
%     while sum_diff > tolerance
    delta = lambda / sqrt(1+lambda^2);
    Delta = sigma * delta;

    [V, W] = trunc_norm_moments(delta / sigma  * (s-u), sqrt(1-delta^2));
    sum_Rs = sum(Rs);

    u_new  = (sum( Rs .* (s - V  * Delta ) )) / (sum_Rs);
    Delta_new  = (sum( Rs .* V  .* (s - u_new ) ) ...
                 ) / (sum(Rs  .* W ));
    Gamma_new  = ( sum( Rs .* ((s - u_new).^2 ...
                               - 2 * V .* (s - u_new) * Delta_new ...
                               + W * Delta_new^2) ) ...
                 ) / (sum_Rs);
    lambda_new = sign(Delta_new) * sqrt(Delta_new^2 / Gamma_new);
    sigma_new = sqrt(Gamma_new + Delta_new^2);

    sum_diff = u_new + sigma_new + lambda_new - u - sigma - lambda;
    u = u_new;
    sigma = sigma_new;
    lambda = lambda_new;
%     end

end