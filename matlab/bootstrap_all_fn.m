function bootstrap_all_fn(bootn, bootstart, bootratio)

list_species = {
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
% 
% 'HeLa01ng_2'
% 'HeLa1ng'
% 'HeLa10ng'
% 'HeLa50ng'
% 'HeLa100ng'

'c_elegans'
'drosophila'
'e_coli'
'human'
'mouse'

}


n = size(list_species, 1)
for i = 1:n
    species = list_species{i}
%     run parse_psm.m
    load(['test_search/matdata/',species,'_data.mat'])
    bootstraping_fn(species, bootstart, bootn, bootratio)
end

end


