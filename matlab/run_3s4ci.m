method='_3s4ci'
[alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3] = EM3_5(omat',1,-1,-1,-1)
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
theta.theta_i2 = pack_skntheta(u_i2, sigma_i2, lambda_i2);
theta.theta_i3 = pack_skntheta(u_i3, sigma_i3, lambda_i3);

save(paramfile, 'theta');
% save(paramfile, ['alpha', 'beta', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i', 'u_i2', 'sigma_i2', 'lambda_i2', 'u_i3', 'sigma_i3', 'lambda_i3'])

plot_fdr
