minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;


s1 = omat(:,1);
s1 = s1(s1~=0);
s2 = omat(:,2);
s2 = s2(s2~=0);
s3 = omat(:,3);
s3 = s3(s3~=0);

% figure('Position', [10,10,1920,1080]);
figure;
hold on;

nsp = 1;
if (strcmp(method(1:5), '_3s4c'))
    nsp = 2;
end
% axis;
yc = normpdf(x_values, u_c, sigma_c);
yi1 = gumbel_pdf(x_values, u_i, sigma_i);
y1 = alpha*yc + (1-alpha)*yi1;
% y4 = skew_norm_pdf
% subplot(nsp,2,1)
% cla;
% hold on;
plot(x_values,yc*alpha,'LineWidth',2);
plot(x_values,yi1*(1-alpha),'LineWidth',2);
plot(x_values,y1,'LineWidth',2);

histogram(s1,200,'Normalization','pdf');

legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

species_folder = ['test_search/est_results/',species];
if ~ exist(species_folder)
    mkdir(species_folder)
    mkdir([species_folder,'/distplot'])
    mkdir([species_folder,'/fitcurv'])
    mkdir([species_folder,'/fdrcurv'])
end
saveas(gcf,[species_folder,'/distplot/',method,'.png'])

% figure('Position', [10,10,1920,1080]);
figure;
hold on;

% subplot(nsp, 2, 1);
% hold on;

S1 = sort(omat(:,1));
m = size(S1,1);
h1emp = (1:m) / m;
pemp = plot(S1,h1emp);

hc = skew_norm_cdf(x_values, u_c, sigma_c, lambda_c);
hi1 = skew_norm_cdf(x_values, u_i, sigma_i, lambda_i);
h1 = alpha*hc + (1-alpha)*hi1;
pest = plot(x_values, h1);

seq = flip(unique(omat(:,1)));

n = size(seq,1);
for i = 1:n
    s = seq(i);
    fdr = fdr_x(alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
%     fdr
    if fdr > 0.1
        break
    end
end
xline(s);
text(s,0.05,'\leftarrow 1%fdr')
xlabel('-log(EValue)');
ylabel('CDF');
ylim([0,1]);

legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});

% subplot(nsp, 2, 2);
% hold on;
% 
% S2 = sort(omat(:,2));
% m = size(S2,1);
% h2emp = (1:m) / m;
% pemp = plot(S2,h2emp);
% 
% hi2 = skew_norm_cdf(x_values, u_i2, sigma_i2, lambda_i2);
% h2 = alpha*hi1 + beta*hc + (1-alpha-beta)*hi2;
% pest = plot(x_values, h2);
% 
% legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});
% ylim([0,1]);

% if (strcmp(method(1:5), '_3s4c'))
%     subplot(nsp,2,3)
%     hold on;
%     
%     S3 = sort(omat(:,3));
%     m = size(S3,1);
%     h3emp = (1:m) / m;
%     pemp = plot(S3,h3emp);
% 
%     hi3 = skew_norm_cdf(x_values, u_i3, sigma_i3, lambda_i3);
%     h3 = (alpha+beta)*hi2 + (1-alpha-beta)*hi3;
%     pest = plot(x_values, h3);
% 
%     legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});
%     ylim([0,1]);
% end


saveas(gcf,[species_folder,'/fitcurv/',method,'.png'])
