
method = '_1s2ca';

load([results_dir,species,'/params/',method,'.mat'])

alpha = theta.alpha;
u_c = theta.theta_c.u;
sigma_c = theta.theta_c.sigma;
a_i = theta.theta_i.a;
b_i = theta.theta_i.b;
gamma_i = theta.theta_i.gamma;


llfile = [ll_folder,method,'.mat'];

ll_1 = func_ll2_1a(s1, alpha, u_c, sigma_c, a_i, b_i, gamma_i) 

save(llfile,'ll_1');

calc_mean_gamma