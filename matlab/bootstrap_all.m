clear

list_species = {
% 'A.thaliana'
% % 'C.elegans'
% 'D.melanogaster'
'E.coli'
% 'H.sapiens2'
% 'H.sapiens3'
% 'M.musculus'
% 'M.musculus2'
% 'M.musculus3'
% % 'S.cerevisiae'
% 'S.cerevisiae2'
% 'S.cerevisiae3'
}


n = size(list_species, 1)
for i = 1:n
    species = list_species{i}
%     run parse_psm.m
    load(['test_search/matdata/',species,'_data.mat'])
    bootstraping
end



