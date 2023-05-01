method='_1s2cg'
[alpha, u_c, sigma_c, u_i, sigma_i] = EM2_1g(omat')
plot_dist_gumbel

species_folder = [results_folder,species]
param_folder = [species_folder,'/params/'];
if ~exist(param_folder)
    mkdir(param_folder)
end
paramfile = [param_folder,method,'.mat'];

theta.alpha = alpha;
theta.theta_c = pack_theta(u_c, sigma_c);
theta.theta_i = pack_theta(u_i, sigma_i);

save(paramfile, 'theta');
% save(paramfile, ['alpha', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i'])

plot_fdr
