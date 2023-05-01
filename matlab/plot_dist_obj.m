
% results_folder = 'test_search/est_results/';

species_folder = [results_folder,species];

figure
hold on
plot(obj.s1, obj.pdfC,'LineWidth',2)
plot(obj.s1, obj.pdfI1,'LineWidth',2)
plot(obj.s1, obj.pdfM,'LineWidth',2)
histogram(obj.s1,200, 'Normalization', 'pdf')
legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});
xlabel('-log(E)')

saveas(gcf,[species_folder,'/distplot/',method,'.png'])


figure
hold on
pest = plot(obj.s1, obj.cdf);
pemp = plot(obj.s1, obj.cdfTrue);

xline(obj.t1p);
text(obj.t1p,0.05,'\leftarrow 1%fdr')
xlabel('-log(EValue)');
ylabel('CDF');
ylim([0,1]);

legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});
xlabel('-log(E)')

saveas(gcf,[species_folder,'/fitcurv/',method,'.png'])

sdcdf = obj.deltaCdf;
sddir = [species_folder, '/sdcdf/'];

save([sddir, method, '.mat'], 'sdcdf');

ll_folder = [species_folder,'/ll/'];
if ~exist(ll_folder)
    mkdir(ll_folder)
end

llfile = [ll_folder,method,'.mat'];
ll_1 = obj.ll;
save(llfile,'ll_1');


