clear

list_species = {
% 'A.thaliana'
% % 'C.elegans'
% 'D.melanogaster'
% 'E.coli'
'H.sapiens2'
% 'H.sapiens3'
% 'M.musculus'
% 'M.musculus2'
% 'M.musculus3'
% % 'S.cerevisiae'
% 'S.cerevisiae2'
% 'S.cerevisiae3'
}
data_dir = 'test_search/matdata_pride/';
results_folder = 'test_search/est_results_pride/';


% species

% results_folder = 'test_search/fragger_results/';

list_methods = {
    '_1s2ca';
    '_1s2c';
    '_2s3ci';
%     '_2s3ct';
};

% list_species = {
%     'c_elegans'
%     'drosophila'
%     'e_coli'
%     'human'
%     'mouse'
%     'yeast'
% %     'human_hcd',
% %     'mouse_hcd',
% };
% data_dir = 'test_search/matdata_nist/';
% results_folder = 'test_search/est_results_nist/';
% xlim_high = 65;
% plotcdf = true;
% plotthres = false;

% list_species = {
%     'HeLa01ng'
%     'HeLa1ng'
%     'HeLa10ng'
%     'HeLa50ng'
%     'HeLa100ng'
% };
% data_dir = 'test_search/matdata_hela/';
% results_folder = 'test_search/est_results_hela/';
% xlim_high = 45;
% plotcdf = true;
% plotthres = true;

figw = 240;
figh = 120;

c_color = [0.8500 0.3250 0.0980];
i1_color = [0.9290 0.6940 0.1250];
i2_color = [0 0.4470 0.7410];
mix_color = [0.4940 0.1840 0.5560];

cdfcolor = [0.6350 0.0780 0.1840];
legending = false;
legending2 = false;

linewidth = 1;

n = size(list_species, 1);

n_m = size(list_methods, 1);

for i = 1:n
    species = list_species{i};

    for me_j = 1:n_m
%         run parse_psm.m
        method = list_methods{me_j};

        species_folder = [results_folder,species];

        load([data_dir, species, '_data.mat'])
        load([species_folder, '/params/', method, '.mat'])

        if strcmp(method, '_1s2c')
            alpha = theta.alpha;
            u_c = theta.theta_c.u;
            sigma_c = theta.theta_c.sigma;
            lambda_c = theta.theta_c.lambda;
            u_i = theta.theta_i.u;
            sigma_i = theta.theta_i.sigma;
            lambda_i = theta.theta_i.lambda;
        elseif strcmp(method, '_1s2ca')
            alpha = theta.alpha;
            u_c = theta.theta_c.u;
            sigma_c = theta.theta_c.sigma;
            a_i = theta.theta_i.a;
            b_i = theta.theta_i.b;
            gamma_i = theta.theta_i.gamma;
        elseif strcmp(method, '_2s3ci') || strcmp(method, '_2s3ct')
            alpha = theta.alpha;
            beta = theta.beta;
            [u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
            [u_i, sigma_i, lambda_i] = unpack_skntheta(theta.theta_i);
            [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta.theta_i2);
        end

        if strcmp(method, '_1s2ca')
            plot_dist_gamma_old
        else
            plot_dist_old
        end

    end
    close all
end


