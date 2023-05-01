
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear
% 
% list_species = {
% 'A.thaliana'
% % 'C.elegans'
% 'D.melanogaster'
% 'E.coli'
% 'H.sapiens2'
% 'H.sapiens3'
% 'M.musculus'
% 'M.musculus2'
% 'M.musculus3'
% % 'S.cerevisiae'
% 'S.cerevisiae2'
% 'S.cerevisiae3'
% };

% list_species = {
%     'HeLa01ng'
%     'HeLa1ng'
%     'HeLa10ng'
%     'HeLa50ng'
%     'HeLa100ng'
%     'HeLa01ng.2'
%     'HeLa1ng.2'
%     'HeLa10ng.2'
%     'HeLa50ng.2'
%     'HeLa100ng.2'
%     'HeLa01ng.3'
%     'HeLa1ng.3'
%     'HeLa10ng.3'
%     'HeLa50ng.3'
%     'HeLa100ng.3'
% };
list_species = {
%     'c_elegans'
%     'drosophila'
%     'e_coli'
%     'human'
%     'mouse'
    'yeast'
};
% species = 'M.musculus'
% species = 'H.sapiens'
% species = 'H.sapiens2'
% species = 'H.sapiens3'
% species = 'H.sapiens4'
% species = 'C.elegans'
% species = 'D.melanogaster'
% species = 'S.cerevisiae'
% species = 'S.cerevisiae2'
% species = 'S.cerevisiae3'
% species = 'E.coli'
% species = 'A.thaliana'

% species

list_methods = {
    '_1s2ca';
    '_1s2c';
    '_2s3ci';
%     '_2s3ct';
};

data_dir = 'test_search/matdata_nist/';

% results_dir = 'test_search/est_results/';
results_dir = 'test_search/est_results_nist/';

n = size(list_species, 1);

m = size(list_methods, 1);

for i = 1:n
    species = list_species{i};

    species_folder = [results_dir,species];
    param_folder = [species_folder,'/params/'];
    ll_folder = [species_folder,'/ll/'];
    if ~exist(ll_folder)
        mkdir(ll_folder)
    end
    mean_folder = [species_folder,'/mean/'];
    if ~exist(mean_folder)
        mkdir(mean_folder)
    end
    for j = 1:m
        method = list_methods{j};

        load([data_dir,species,'_data.mat'])

        S = omat';
        s1 = S(1,:);
        s1 = s1(s1~=0);
        s2 = S(2,:);
        s2 = s2(s2~=0);
        
        run(['ll',method])
    end
end
