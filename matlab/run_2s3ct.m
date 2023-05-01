method='_2s3ct'
[alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2] = EM3_7(omat',1,1)
plot_dist

species_folder = [results_folder,species];

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
theta.beta = beta;
theta.theta_c = pack_skntheta(u_c, sigma_c, lambda_c);
theta.theta_i = pack_skntheta(u_i, sigma_i, lambda_i);
theta.theta_i2 = pack_skntheta(u_i2, sigma_i2, lambda_i2);

save(paramfile, 'theta');
% save(paramfile, ['alpha', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i'])

plot_fdr
