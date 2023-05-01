
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear
% 
% list_species = {
% % 'A.thaliana'
% % % 'C.elegans'
% % 'D.melanogaster'
% % 'E.coli'
% 'H.sapiens2'
% % 'H.sapiens3'
% % 'M.musculus'
% % 'M.musculus2'
% % 'M.musculus3'
% % % 'S.cerevisiae'
% % 'S.cerevisiae2'
% % 'S.cerevisiae3'
% }
% data_dir = 'test_search/matdata_pride/';
% 
% results_folder = 'test_search/est_results_pride/';
% xlim_low = 0;
% xlim_high = 65;
% bin_width = 0.5;
% plotcdf = true;
% plotthres = false;

% list_species = {
% 'HeLa01ng'
% 'HeLa1ng'
% 'HeLa10ng'
% 'HeLa50ng'
% 'HeLa100ng'
% 'HeLa01ng.2'
% 'HeLa1ng.2'
% 'HeLa10ng.2'
% 'HeLa50ng.2'
% 'HeLa100ng.2'
% 'HeLa01ng.3'
% 'HeLa1ng.3'
% 'HeLa10ng.3'
% 'HeLa50ng.3'
% 'HeLa100ng.3'
% }
% data_dir = 'test_search/matdata_hela/';
% 
% results_folder = 'test_search/est_results_hela';
% xlim_low = 0;
% xlim_high = 65;
% bin_width = 0.5;
% plotcdf = true;
% plotthres = false;


% list_species = {
% %     'c_elegans'
% %     'drosophila'
% %     'e_coli'
% %     'human'
% %     'mouse'
%     'yeast'
% }
% 
% data_dir = 'test_search/matdata_nist/';
% 
% results_folder = 'test_search/est_results_nist';
% xlim_low = 0;
% xlim_high = 65;
% bin_width = 0.5;
% plotcdf = true;
% plotthres = false;

list_species = {
    'ecoli'
}
data_dir = 'test_search/matdata_xl/';

results_folder = 'test_search/est_results_xl/';
xlim_low = -8;
xlim_high = 6;
bin_width = 0.5;
plotcdf = true;
plotthres = true;


figw = 240;
figh = 120;

c_color = [0.8500 0.3250 0.0980];
i1_color = [0.9290 0.6940 0.1250];
i2_color = [0 0.4470 0.7410];
ic_color = [0 0.4470 0.7410];
mix_color = [0.4940 0.1840 0.5560];

cdfcolor = [0.6350 0.0780 0.1840];
legending = true;
legending2 = false;

linewidth = 1;

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
%     run parse_psm.m
    load([data_dir,species,'_data.mat'])
    run_all
end
