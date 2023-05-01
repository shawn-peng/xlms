function condition_satisfied = check_skewnorm_pdf_higher(x, ...
    u_1, sigma_1, lambda_1, ...
    u_2, sigma_2, lambda_2)

p1 = skew_norm_pdf(x, u_1, sigma_1, lambda_1);
p2 = skew_norm_pdf(x, u_2, sigma_2, lambda_2);

condition_satisfied = all(p1 >= p2);

end