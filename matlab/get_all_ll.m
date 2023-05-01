list_species = {
'M.musculus'
'M.musculus2'
'M.musculus3'
'H.sapiens2'
'H.sapiens3'
% 'C.elegans'
'D.melanogaster'
% 'S.cerevisiae'
'S.cerevisiae2'
'S.cerevisiae3'
'E.coli'
'A.thaliana'
}

method = '_2s3c'

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
    load(['test_search/matdata/',species,'_data.mat'])
    load(['test_search/est_results/',species,'/params/',method,'.mat'])
    alpha = theta.alpha;
    beta = theta.beta;
    [u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
    [u_i1, sigma_i1, lambda_i1] = unpack_skntheta(theta.theta_i);
    [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta.theta_i2);
%     [u_i3, sigma_i3, lambda_i3] = unpack_skntheta(theta.theta_i3);
    S = omat';
    s1 = S(1,:);
    s1 = s1(s1~=0);
    s2 = S(2,:);
    s2 = s2(s2~=0);
    s3 = S(3,:);
    s3 = s3(s3~=0);
%     ll = func_ll2_1(s1, alpha, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1);
%     ll = func_ll3_5(s1, s2, s3, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    ll = func_ll3(s1, s2, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2);
    ll = ll / size(s1,2);
    disp({species; ll});
end

