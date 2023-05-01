method='_1s2c'
species_folder = [results_folder,species]
% [alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = EM2_1(omat',1,1)
em_1s2c_with_4init
plot_dist

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

theta.alpha = alpha;
theta.theta_c = pack_skntheta(u_c, sigma_c, lambda_c);
theta.theta_i = pack_skntheta(u_i, sigma_i, lambda_i);

save(paramfile, 'theta');
% save(paramfile, ['alpha', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i'])

plot_fdr
