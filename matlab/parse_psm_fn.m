function parse_psm_fn(psm_dir, species)

species

filename = [psm_dir, species, '_nod.tsv']
x = tdfread(filename);


% evalue_field = 'EValue'
% evalue_field = 'SpecEValue'

n = size(x.SpecID, 1);

ns = 2;

curr_scannum = x.ScanNum(1);
sa = [];
smat = [];
mat = [];
slen = zeros(max(x.ScanNum),1);
% slen = [];
for i = 1:n
%     scannum = x.ScanNum(i);
%     if scannum == -1
%         print('scannum == -1')
        specid = x.SpecID(i,:);
        kvpairs = split(specid);
        scannum = split(kvpairs(1), '=');
        scannum = str2num(scannum{2})+1;
%     end
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
    if slen(scannum) > 0 && x.SpecEValue(i) == smat(scannum, slen(scannum))
        continue
    end
    slen(scannum) = slen(scannum) + 1;
%     mat(scannum, slen(scannum)) = x.MSGFScore(i);
    smat(scannum, slen(scannum)) = x.SpecEValue(i);
    mat(scannum, slen(scannum)) = -log(x.SpecEValue(i));
end

specind = find(slen>=3);
mat3 = mat(specind,1:3);
smat = smat(slen>0,:);
omat = mat(slen>0,:);
% emptyind = find(slen==1);
% mat(emptyind,2) = -inf;
specind = find(slen>=2);
% mat = mat(specind,:);
mat2 = mat(specind,1:2);


save(['test_search/matdata/',species,'_data.mat'], 'species', 'omat', 'smat', 'mat', 'mat2', 'mat3', 'slen');
% save(['test_search/matdata_openms/',species,'_data.mat'], 'species', 'omat', 'smat', 'mat', 'mat2', 'mat3', 'slen');

