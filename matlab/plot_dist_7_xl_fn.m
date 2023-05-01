function plot_dist_7_xl_fn(fig, mat, method, ws, ...
u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i, sigma_i, lambda_i, ...
u_c2, sigma_c2, lambda_c2, u_ic2, sigma_ic2, lambda_ic2, u_i12, sigma_i12, lambda_i12, u_i22, sigma_i22, lambda_i22)
minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;

mat = mat';
omat = mat;
% mat3 = mat;
mat2 = mat;

s1 = mat(:,1);
% s1 = s1(s1~=0);
s2 = mat(:,2);
% s2 = s2(s2~=0);
% s3 = mat(:,3);
% s3 = s3(s3~=0);

wc = num2cell(ws);
[w1, w2, w3, w4, w5, w6, w7] = wc{:};

w1c = (w1+w2);
w1ic = (w3+w5);
w1i1 = (w4+w6+w7);
w2c = (w3+w4);
w2ic = (w1+w6);
w2i1 = (w2+w5);
w2i2 = w7;

% axis;

yc = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
yic = skew_norm_pdf(x_values, u_ic, sigma_ic, lambda_ic);
yi1 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);

% y1 = alpha*yc + beta*yi1 + (1-alpha-beta)*yic;
y1 = (w1+w2)*yc + (w3+w5)*yic + (w4+w6+w7)*yi1;
% y4 = skew_norm_pdf

figure(fig);
ax_1 = subplot(2,1,1);
cla;
hold on;
plot(ax_1,x_values,yc*w1c,'LineWidth',2);
plot(ax_1,x_values,yic*w1ic,'LineWidth',2, 'LineStyle', '--');
plot(ax_1,x_values,yi1*w1i1,'LineWidth',2);
plot(ax_1,x_values,y1,'LineWidth',2);

histogram(ax_1,s1,100,'Normalization','pdf', 'FaceColor', 'none');

legend(ax_1,{'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

xlim([minS, maxS]);

figure(fig);
ax_2 = subplot(2,1,2);
cla;
hold on;

yc2 = skew_norm_pdf(x_values, u_c2, sigma_c2, lambda_c2);
yic2 = skew_norm_pdf(x_values, u_ic2, sigma_ic2, lambda_ic2);
yi12 = skew_norm_pdf(x_values, u_i12, sigma_i12, lambda_i12);
yi22 = skew_norm_pdf(x_values, u_i22, sigma_i22, lambda_i22);
y2 = (w1+w6)*yic2 + (w2+w5)*yi12 + (w3+w4)*yc2 + w7*yi22;
plot(ax_2,x_values,yc2*w2c,'LineWidth',2);
plot(ax_2,x_values,yic2*w2ic,'LineWidth',2, 'LineStyle', '--');
plot(ax_2,x_values,yi12*w2i1,'LineWidth',2);
plot(ax_2,x_values,yi22*w2i2,'LineWidth',2, 'LineStyle', ':');
plot(ax_2,x_values,y2,'LineWidth',2);
histogram(ax_2,s2,100,'Normalization','pdf', 'FaceColor', 'none');
legend(ax_2,{'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'dist\_i2'; 'mixture'; 'hist\_second'});

xlim([minS, maxS]);

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


