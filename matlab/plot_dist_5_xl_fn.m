function plot_dist_5_xl_fn(mat, method, ws, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_i, sigma_i, lambda_i, ...
    u_ic2, sigma_ic2, lambda_ic2, ...
    u_i2, sigma_i2, lambda_i2)

minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;

mat = mat';
omat = mat;
% mat3 = mat;
mat2 = mat;

s1 = mat(:,1);
s1 = s1(s1~=0);
s2 = mat(:,2);
s2 = s2(s2~=0);
% s3 = mat(:,3);
% s3 = s3(s3~=0);

wc = num2cell(ws);
[w1, w2, w3, w4, w5, w6, w7, w8] = wc{:};

w1c = (w1+w2);
w1ic = (w3+w5+w7);
w1i1 = (w4+w6+w8);
w2c = (w3+w4);
w2ic = (w1+w6);
w2i1 = (w2+w5);
w2ic2 = w7;
w2i2 = w8;

% axis;

yc = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
yic = skew_norm_pdf(x_values, u_ic, sigma_ic, lambda_ic);
yi1 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);
yic2 = skew_norm_pdf(x_values, u_ic2, sigma_ic2, lambda_ic2);
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);

% y1 = alpha*yc + beta*yi1 + (1-alpha-beta)*yic;
y1 = w1c*yc + w1ic*yic + w1i1*yi1;
% y4 = skew_norm_pdf

subplot(2,1,1)
cla;
hold on;
plot(x_values,yc*w1c,'LineWidth',2);
plot(x_values,yic*w1ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi1*w1i1,'LineWidth',2);
plot(x_values,y1,'LineWidth',2);

xlim([minS, maxS]);
histogram(s1,100,'Normalization','pdf', 'FaceColor', 'none');

legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

subplot(2,1,2)
cla;
hold on;
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);
y2 = w2ic*yic + w2i1*yi1 + w2c*yc + w2ic2*yic2 + w2i2*yi2;
plot(x_values,yc*w2c,'LineWidth',2);
plot(x_values,yic*w2ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi1*w2i1,'LineWidth',2);
plot(x_values,yic2*w2ic2,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi2*w2i2,'LineWidth',2, 'LineStyle', ':');
plot(x_values,y2,'LineWidth',2);
histogram(s2,100,'Normalization','pdf', 'FaceColor', 'none');
xlim([minS, maxS]);
legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'dist\_half\_incorrect2'; 'dist\_incorrect2';  'mixture'; 'hist\_second'});

hold off;
subplot(2,1,1)

end


