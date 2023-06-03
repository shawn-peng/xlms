function condition_satisfied = check_skewnorm_cdf_higher(x, ...
    u_1, sigma_1, lambda_1, ...
    u_2, sigma_2, lambda_2, ...
    fuzzy)

p1 = skew_norm_cdf(x, u_1, sigma_1, lambda_1);
p2 = skew_norm_cdf(x, u_2, sigma_2, lambda_2);

condition_satisfied = all((p2 - p1) >= -fuzzy);

end