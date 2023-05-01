function bootstraping_fn(species, bootstart, nboots, boot_ratio)
% bootstart = 181;
% nboots = 20;
% boot_ratio = .1;
bootind_dir = 'test_search/est_results/bootstrap_index/';
if ~exist(bootind_dir)
	mkdir(bootind_dir)
end
for bootstrap_num = bootstart:bootstart+nboots-1
    loaddata
    n = size(mat,1);
    boot_n = round(n * boot_ratio);
    ind = 1:n;
    boot_ind = datasample(ind, boot_n);
	writematrix(boot_ind, [bootind_dir, num2str(bootstrap_num), '.csv'], 'Delimiter', 'tab');

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

end
