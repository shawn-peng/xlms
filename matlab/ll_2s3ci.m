
method = '_2s3ci';

load([results_dir,species,'/params/',method,'.mat'])

alpha = theta.alpha;
beta = theta.beta;
[u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
[u_i, sigma_i, lambda_i] = unpack_skntheta(theta.theta_i);
[u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta.theta_i2);


llfile = [ll_folder,method,'.mat'];

ll_1 = func_ll2_1(s1, alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i)

ll_12 = func_ll3(s1, s2, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2)

save(llfile,'ll_1','ll_12');

calc_mean

