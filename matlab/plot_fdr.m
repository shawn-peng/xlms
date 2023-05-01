

curv = [0,0,0];

seq = flip(unique(omat(:,1)));

n = size(seq,1);
for i = 1:n
    s = seq(i);
    fdr = fdr_x(alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
%     fdr
    if fdr > 0.1
        break
    end
    curv(i,:) = [fdr, sum(omat(:,1)>s), s];
end

figure();
plot(curv(:,1), curv(:,2));

% csvwrite('twomix.csv', curv);
species_folder = [results_folder,species];
fdr_folder = [species_folder,'/fdr/'];
if exist('bootstrap_num')
    fdr_folder = [fdr_folder,method,'/bootstrap/'];
end
if ~exist(fdr_folder)
    mkdir(fdr_folder)
end
if exist('bootstrap_num')
    fdrfile = [fdr_folder,num2str(bootstrap_num),'.csv'];
else
    fdrfile = [fdr_folder,method,'.csv'];
end

dlmwrite(fdrfile, curv, 'delimiter', ',', 'precision', '%25.20f');
