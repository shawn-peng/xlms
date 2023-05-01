
% x = tdfread('test_search/Adult_Adrenalgland_Gel_Elite_49_f01_nod.tsv');


% species = 'M.musculus'
% species = 'M.musculus2'
% species = 'M.musculus3'
species = 'H.sapiens2'
% species = 'H.sapiens3'
% species = 'C.elegans'
% species = 'D.melanogaster'
% species = 'S.cerevisiae'
% species = 'S.cerevisiae2'
% species = 'S.cerevisiae3'
% species = 'E.coli'
% species = 'A.thaliana'

% method = '_3s4c'
% method = '_3mix'
% method = '_2mix'


% x = tdfread('test_search/isb02_t3_nodecoy.tsv');
% x = tdfread('test_search/pride/H.sapiens_nod.tsv');
% x = tdfread('test_search/pride/M.musculus_nod.tsv');
% x = tdfread('test_search/pride/C.elegans_nod.tsv');
% x = tdfread('test_search/pride/D.melanogaster_nod.tsv');
% x = tdfread('test_search/pride/S.cerevisiae_nod.tsv');
% x = tdfread('test_search/pride/P.aeruginosa_nod.tsv');
% x = tdfread('test_search/pride/E.coli_nod.tsv');
% x = tdfread('test_search/pride/A.thaliana_nod.tsv');
species
% filename = sprintf('test_search/pride/%s_nod.tsv', species);
filename = sprintf('test_search/pride/fragger/%s_nod.tsv', species);
x = tdfread(filename);


% evalue_field = 'EValue'
% evalue_field = 'SpecEValue'

n = size(x.scannum, 1);

ns = 2;

curr_scannum = x.scannum(1);
sa = [];
smat = [];
mat = [];
slen = zeros(max(x.scannum),1);
% slen = [];
for i = 1:n
    scannum = x.scannum(i);
    if mod(i,100) == 0
        disp(sprintf('%d/%d',i,n));
    end
    curr_scannum = scannum;
    if scannum > size(slen)
        slen(scannum) = 0;
    end
%     if slen(scannum) >= 2
%         continue
%     end
    if slen(scannum) > 0 && x.expectscore(i) == smat(scannum, slen(scannum))
        continue
    end

    slen(scannum) = slen(scannum) + 1;
%     mat(scannum, slen(scannum)) = x.MSGFScore(i);
    smat(scannum, slen(scannum)) = x.expectscore(i);
    mat(scannum, slen(scannum)) = -log(x.expectscore(i));
end

flags = mat(:,1)==inf;
slen(flags) = 0;
mat(flags,:) = 0;

specind = find(slen>=3);
mat3 = mat(specind,1:3);
smat = smat(specind,1:3);
omat = mat(slen>0,:);
% emptyind = find(slen==1);
% mat(emptyind,2) = -inf;
specind = find(slen>=2);
% mat = mat(specind,:);
mat2 = mat(specind,1:2);


save(['test_search/matdata/fragger/',species,'_data.mat'], 'species', 'omat', 'smat', 'mat', 'mat2', 'mat3', 'slen');

