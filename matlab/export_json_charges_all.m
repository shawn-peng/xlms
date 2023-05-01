
clear

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
% }
list_species = {
    'HeLa01ng'
    'HeLa1ng'
    'HeLa10ng'
    'HeLa50ng'
    'HeLa100ng'
    'HeLa01ng.2'
    'HeLa1ng.2'
    'HeLa10ng.2'
    'HeLa50ng.2'
    'HeLa100ng.2'
    'HeLa01ng.3'
    'HeLa1ng.3'
    'HeLa10ng.3'
    'HeLa50ng.3'
    'HeLa100ng.3'
};
% list_species = {
%     'c_elegans'
%     'drosophila'
%     'e_coli'
%     'human'
%     'mouse'
% }

json_dir = 'est_results/json/';
if ~exist(json_dir)
    mkdir(json_dir)
end

% method = '_1s2c';

results_folder = 'test_search/est_results/';


list_methods = {
    '_1s2ca';
    '_1s2c';
    '_2s3ci';
%     '_2s3ct';
};

n_sp = size(list_species, 1);

n_m = size(list_methods, 1);

results_arr = {};

for sp_i = 1:n_sp
    species = list_species{sp_i};
    export_results_json
    for c = 2:3
        species = list_species{sp_i};
        species = [species, '.c', num2str(c)];
        export_results_json
    end
end


