function plot_dist_10_xl_fn(fig, mat, method, ws, ...
u_c, sigma_c, lambda_c, ...
u_ic, sigma_ic, lambda_ic, ...
u_i, sigma_i, lambda_i, ...
u_c2ic, sigma_c2ic, lambda_c2ic, ...
u_c2i1, sigma_c2i1, lambda_c2i1, ...
u_ic2c, sigma_ic2c, lambda_ic2c, ...
u_ic2i1, sigma_ic2i1, lambda_ic2i1, ...
u_i12c, sigma_i12c, lambda_i12c, ...
u_i12ic, sigma_i12ic, lambda_i12ic, ...
u_i22, sigma_i22, lambda_i22)
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
w1ic = (w3+w4);
w1i1 = (w5+w6+w7);
w2ic2c = w1;
w2i12c = w2;
w2c2ic = w3;
w2i12ic = w4;
w2c2i1 = w5;
w2ic2i1 = w6;
w2i22 = w7;

% axis;

yc = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
yic = skew_norm_pdf(x_values, u_ic, sigma_ic, lambda_ic);
yi1 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);

% y1 = alpha*yc + beta*yi1 + (1-alpha-beta)*yic;
y1 = (w1+w2)*yc + (w3+w4)*yic + (w5+w6+w7)*yi1;
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

yc2ic = skew_norm_pdf(x_values, u_c2ic, sigma_c2ic, lambda_c2ic);
yc2i1 = skew_norm_pdf(x_values, u_c2i1, sigma_c2i1, lambda_c2i1);
yic2c = skew_norm_pdf(x_values, u_ic2c, sigma_ic2c, lambda_ic2c);
yic2i1 = skew_norm_pdf(x_values, u_ic2i1, sigma_ic2i1, lambda_ic2i1);
yi12c = skew_norm_pdf(x_values, u_i12c, sigma_i12c, lambda_i12c);
yi12ic = skew_norm_pdf(x_values, u_i12ic, sigma_i12ic, lambda_i12ic);
yi22 = skew_norm_pdf(x_values, u_i22, sigma_i22, lambda_i22);

y2 = w2c2ic*yc2ic + w2c2i1*yc2i1 + w2ic2c*yic2c + w2ic2i1*yic2i1 + w2i12c*yi12c + w2i12ic*yi12ic + w2i22*yi22;

plot(ax_2,x_values,yc2i1*w2c2i1,'LineWidth',2);
plot(ax_2,x_values,yic2c*w2ic2c + yic2i1*w2ic2i1 + yc2ic*w2c2ic,'LineWidth',2, 'LineStyle', '--');
plot(ax_2,x_values,yi12c*w2i12c + yi12ic*w2i12ic,'LineWidth',2);
plot(ax_2,x_values,yi22*w2i22,'LineWidth',2, 'LineStyle', ':');
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


