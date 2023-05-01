method='_1s2ca'
[alpha, u_c, sigma_c, a_i, b_i, gamma_i] = EM2_1a(omat')
plot_dist_gamma

species_folder = [results_folder,species]
species_folder
if ~exist('bootstrap_num')
    param_folder = [species_folder,'/params/'];
    paramfile = [param_folder,method,'.mat'];
else
    param_folder = [species_folder,'/params/',method,'/bootstrap/'];
    paramfile = [param_folder,num2str(bootstrap_num),'.mat'];
end
if ~exist(param_folder)
    mkdir(param_folder)
end

theta = {};
theta.alpha = alpha;
theta.theta_c.u = u_c;
theta.theta_c.sigma = sigma_c;
theta.theta_i.a = a_i;
theta.theta_i.b = b_i;
theta.theta_i.gamma = gamma_i;


save(paramfile, 'theta');
% save(paramfile, ['alpha', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i'])

plot_fdr_a
