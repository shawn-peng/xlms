
load([results_dir,species,'/params/',method,'.mat'])

u_c = theta.theta_c.u;
sigma_c = theta.theta_c.sigma;
a_i = theta.theta_i.a;
b_i = theta.theta_i.b;
gamma_i = theta.theta_i.gamma;


% p_c = skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
% p_i = skew_norm_pdf(s1, u_i, sigma_i, lambda_i);

E_c = u_c
E_i = a_i / b_i + gamma_i

mean_file = [mean_folder,method,'.mat'];

save(mean_file, 'E_c', 'E_i')
