
load([results_dir,species,'/params/',method,'.mat'])

[u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
[u_i, sigma_i, lambda_i] = unpack_skntheta(theta.theta_i);


% p_c = skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
% p_i = skew_norm_pdf(s1, u_i, sigma_i, lambda_i);

E_c = u_c + sigma_c * lambda_c / sqrt(1 + lambda_c^2) * sqrt(2/pi)
E_i = u_i + sigma_i * lambda_i / sqrt(1 + lambda_i^2) * sqrt(2/pi)

mean_file = [mean_folder,method,'.mat'];

save(mean_file, 'E_c', 'E_i')
