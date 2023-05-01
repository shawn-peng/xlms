
species = 'A.thaliana'
% species = 'C.elegans'
% species = 'D.melanogaster'
% species = 'E.coli'
% species = 'H.sapiens'
% species = 'H.sapiens2'
% species = 'H.sapiens3'
% species = 'M.musculus'
% species = 'M.musculus2'
% species = 'M.musculus3'
% species = 'S.cerevisiae'
% species = 'S.cerevisiae2'
% species = 'S.cerevisiae3'

method = '_2s3ci'

param_file = ['test_search/est_results/',species,'/params/',method,'.mat'];
load(param_file)

[u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
[u_i, sigma_i, lambda_i] = unpack_skntheta(theta.theta_i);
if ~strcmp(method(1:5), '_1s2c')
    [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta.theta_i2);
end
if strcmp(method(1:5), '_3s4c')
    [u_i3, sigma_i3, lambda_i3] = unpack_skntheta(theta.theta_i3);
end

[u_c, u_i, u_i2, u_i3]

