list_species = {
% 'M.musculus'
% 'M.musculus2'
% 'M.musculus3'
% 'H.sapiens2'
% 'H.sapiens3'
% 'C.elegans'
% 'D.melanogaster'
% 'S.cerevisiae'
% 'S.cerevisiae2'
% 'S.cerevisiae3'
% 'E.coli'
% 'A.thaliana'

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

% 'c_elegans'
% 'drosophila'
% 'e_coli'
% 'human'
% 'mouse'
}

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
    run parse_psm.m
%     run parse_nist_psm.m
end
