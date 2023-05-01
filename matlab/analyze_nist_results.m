
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear

list_species = {
% 'c_elegans'
% 'drosophila'
% 'e_coli'
% 'human'
'mouse'
}

% species

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
%     run parse_psm.m
    load(['test_search/matdata/nist/',species,'_data.mat'])
    run_all
end
