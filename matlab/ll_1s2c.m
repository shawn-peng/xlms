
method = '_1s2c';

load([results_dir,species,'/params/',method,'.mat'])

alpha = theta.alpha;
u_c = theta.theta_c.u;
sigma_c = theta.theta_c.sigma;
lambda_c = theta.theta_c.lambda;
u_i = theta.theta_i.u;
sigma_i = theta.theta_i.sigma;
lambda_i = theta.theta_i.lambda;


llfile = [ll_folder,method,'.mat'];

ll_1 = func_ll2_1(s1, alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i)

save(llfile,'ll_1');

calc_mean
