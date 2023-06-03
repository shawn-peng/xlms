function plot_dist_i_xl_ic2_fn(mat, method, ws, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_ic2, sigma_ic2, lambda_ic2, ...
    u_i, sigma_i, lambda_i, ...
    u_i2, sigma_i2, lambda_i2)
minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
% x_values = minS:0.1:maxS;
x_values = 0:0.1:maxS;

mat = mat';
omat = mat;
% mat3 = mat;
mat2 = mat;

s1 = mat(:,1);
s1 = s1(s1~=0);
% s1 = s1(s1>0);
s2 = mat(:,2);
s2 = s2(s2~=0);
% s2 = s2(s2>0);
% s3 = mat(:,3);
% s3 = s3(s3~=0);

wc = num2cell(ws);
[w1c, w1ic, w1i1, w2c, w2ic, w2ic2, w2i1, w2i2] = wc{:};

% axis;
yc  = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
yic = skew_norm_pdf(x_values, u_ic, sigma_ic, lambda_ic);
yic2 = skew_norm_pdf(x_values, u_ic2, sigma_ic2, lambda_ic2);
yi1 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);
% y1 = alpha*yc + beta*yi1 + (1-alpha-beta)*yic;
y1  = w1c*yc + w1ic*yic + w1i1*yi1;
% y4 = skew_norm_pdf


ax1 = subplot(2,1,1);
cla;
hold on;
plot(x_values,yc*w1c,'LineWidth',2);
plot(x_values,yic*w1ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi1*w1i1,'LineWidth',2);
plot(x_values,y1,'LineWidth',2);

histogram(s1,100,'Normalization','pdf', 'FaceColor', 'none');

xlim([0, maxS]);

legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

% yic_n = normpdf(x_values, u_ic, sigma_ic);
% plot(x_values,yic_n*w1ic,'LineWidth',2, 'LineStyle', '--');
hold off;

ax2 = subplot(2,1,2);
delete(ax2)
ax2 = subplot(2,1,2);
hold on;
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);
y2 = w2ic*yic + w2ic2*yic2 + w2i1*yi1 + w2c*yc + w2i2*yi2;
plot(x_values,yc*w2c,'LineWidth',2);
plot(x_values,yic*w2ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yic2*w2ic2,'LineWidth',2, 'LineStyle', '-.');
plot(x_values,yi1*w2i1,'LineWidth',2);
plot(x_values,yi2*w2i2,'LineWidth',2, 'LineStyle', ':');
plot(x_values,y2,'LineWidth',2);
histogram(s2,100,'Normalization','pdf', 'FaceColor', 'none');
legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_half\_incorrect2'; 'dist\_incorrect'; 'dist\_i2'; 'mixture'; 'hist\_second'});

xlim([0, maxS]);
% plot(x_values,yi2*w2i1,'LineWidth',2, 'LineStyle', ':');

linkaxes([ax1, ax2], 'x');

yyaxis(ax2, 'right');
cla;
% plot_skewnorm_cdf(ax2, x_values, u_i, sigma_i, lambda_i, 1);
% plot_skewnorm_cdf(ax2, x_values, u_i2, sigma_i2, lambda_i2, 1);
hold off;

if (strcmp(method(1:5), '_3s4c'))
    subplot(nsp,2,3)
    cla;
    hold on;
    yi3 = skew_norm_pdf(x_values, u_i3, sigma_i3, lambda_i3);
    y3 = alpha*yi2 + (1-alpha)*yi3;
    plot(x_values,alpha*yi2,'LineWidth',2);
    plot(x_values,(1-alpha)*yi3,'LineWidth',2);
    plot(x_values,y3,'LineWidth',2);
    histogram(s3,100,'Normalization','pdf');
    legend({'dist\_i2'; 'dist\_i3'; 'mixture'; 'hist\_third'});
%     histogram(mat3(:,2),100,'Normalization','pdf');
%     histogram(mat3(:,3),100,'Normalization','pdf');
end

% legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'mixture2'; 'mixture3'; 'hist\_first'; 'hist\_second'});

hold off;

end


