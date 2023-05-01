% clear
% species = 'S.cerevisiae3'
bootstart = 181;
nboots = 20;
boot_ratio = .1;
for bootstrap_num = bootstart:bootstart+nboots-1
    loaddata
    n = size(mat,1);
    boot_n = round(n * boot_ratio);
    ind = 1:n;
    boot_ind = datasample(ind, boot_n);
    mat = mat(boot_ind,:);
%     smat = smat(boot_ind,:);
    slen = slen(boot_ind);

    specind = find(slen>=3);
    mat3 = mat(specind,1:3);
%     smat = smat(specind,1:3);
    omat = mat(slen>0,:);
    % emptyind = find(slen==1);
    % mat(emptyind,2) = -inf;
    specind = find(slen>=2);
    % mat = mat(specind,:);
    mat2 = mat(specind,1:2);
    
    run_all
    close all
end
