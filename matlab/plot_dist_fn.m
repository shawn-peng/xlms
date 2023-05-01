function plot_dist_fn(mat, method, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3)
minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;

mat = mat';
omat = mat;
mat3 = mat;
mat2 = mat;

s1 = mat(:,1);
s1 = s1(s1~=0);
s2 = mat(:,2);
s2 = s2(s2~=0);


nsp = 1;
if (strcmp(method(1:5), '_3s4c'))
    s3 = mat(:,3);
    s3 = s3(s3~=0);
    nsp = 2;
end
% axis;
yc = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
yi1 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);
y1 = alpha*yc + (1-alpha)*yi1;
% y4 = skew_norm_pdf
subplot(nsp,2,1)
cla;
hold on;
plot(x_values,yc*alpha,'LineWidth',2);
plot(x_values,yi1*(1-alpha),'LineWidth',2);
plot(x_values,y1,'LineWidth',2);


histogram(s1,200,'Normalization','pdf', 'FaceColor', 'none');

legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

subplot(nsp,2,2)
cla;
hold on;
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);
y2 = alpha*yi1 + beta*yc + (1-alpha-beta)*yi2;
plot(x_values,yc*beta,'LineWidth',2);
plot(x_values,yi1*alpha,'LineWidth',2);
plot(x_values,y2,'LineWidth',2);
plot(x_values,yi2*(1-alpha-beta),'LineWidth',2);
histogram(s2,100,'Normalization','pdf', 'FaceColor', 'none');
legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'dist\_i2'; 'hist\_second'});

if (strcmp(method(1:5), '_3s4c'))
    subplot(nsp,2,3)
    cla;
    hold on;
    yi3 = skew_norm_pdf(x_values, u_i3, sigma_i3, lambda_i3);
    y3 = alpha*yi2 + (1-alpha)*yi3;
    plot(x_values,alpha*yi2,'LineWidth',2);
    plot(x_values,(1-alpha)*yi3,'LineWidth',2);
    plot(x_values,y3,'LineWidth',2);
    histogram(s3,100,'Normalization','pdf', 'FaceColor', 'none');
    legend({'dist\_i2'; 'dist\_i3'; 'mixture'; 'hist\_third'});
%     histogram(mat3(:,2),100,'Normalization','pdf');
%     histogram(mat3(:,3),100,'Normalization','pdf');
end

% legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'mixture2'; 'mixture3'; 'hist\_first'; 'hist\_second'});

hold off;

% figure;
% hold on;
% 
% S1 = sort(omat(:,1));
% m = size(S1,1);
% h1 = (1:m) / m;
% pemp = plot(S1,h1);
% 
% h11 = skew_norm_cdf(x_values, u_c, sigma_c, lambda_c);
% h12 = skew_norm_cdf(x_values, u_i, sigma_i, lambda_i);
% h1e = alpha*h11 + (1-alpha)*h12;
% pest = plot(x_values, h1e);
% 
% seq = flip(unique(omat(:,1)));
% 
% n = size(seq,1);
% for i = 1:n
%     s = seq(i);
%     fdr = fdr_x(alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
% %     fdr
%     if fdr > 0.1
%         break
%     end
% end
% xline(s);
% text(s,0.05,'\leftarrow 1%fdr')
% xlabel('-log(EValue)');
% ylabel('CDF');
% ylim([0,1]);
% 
% legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});

% saveas(gcf,['test_search/fitting curve/',species,method,'.png'])

end
