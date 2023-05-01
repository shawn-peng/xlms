method='_2s2ci'
[alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = EM2_2(omat',1,1)
plot_dist

species_folder = [results_folder,species]
param_folder = [species_folder,'/params/'];
if ~exist(param_folder)
    mkdir(param_folder)
end
paramfile = [param_folder,method,'.mat'];

theta.alpha = alpha;
theta.beta = beta;
theta.theta_c = pack_skntheta(u_c, sigma_c, lambda_c);
theta.theta_i = pack_skntheta(u_i, sigma_i, lambda_i);

save(paramfile, 'theta');
% save(paramfile, ['alpha', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i'])

plot_fdr
