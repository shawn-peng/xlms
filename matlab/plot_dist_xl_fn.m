function plot_dist_xl_fn(mat, method, ws, ...
u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2)
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
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);

% y1 = alpha*yc + beta*yi1 + (1-alpha-beta)*yic;
y1 = (w1+w2)*yc + (w3+w5)*yic + (w4+w6+w7)*yi1;
% y4 = skew_norm_pdf

subplot(2,1,1)
cla;
hold on;
plot(x_values,yc*w1c,'LineWidth',2);
plot(x_values,yic*w1ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi1*w1i1,'LineWidth',2);
plot(x_values,y1,'LineWidth',2);


histogram(s1,100,'Normalization','pdf', 'FaceColor', 'none');

legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

subplot(2,1,2)
cla;
hold on;
yi2 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);
y2 = (w1+w6)*yic + (w2+w5)*yi1 + (w3+w4)*yc + w7*yi2;
plot(x_values,yc*w2c,'LineWidth',2);
plot(x_values,yic*w2ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi1*w2i1,'LineWidth',2);
plot(x_values,yi2*w2i2,'LineWidth',2, 'LineStyle', ':');
plot(x_values,y2,'LineWidth',2);
histogram(s2,100,'Normalization','pdf', 'FaceColor', 'none');
legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'dist\_i2'; 'mixture'; 'hist\_second'});

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


