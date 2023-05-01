function plot_dist_1i_xl_fn(mat, method, ws, ...
u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
u_i, sigma_i, lambda_i)
minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;

mat = mat';

s1 = mat(:,1);
s1 = s1(s1~=0);
% s3 = mat(:,3);
% s3 = s3(s3~=0);

wc = num2cell(ws);
[w1c, w1ic, w1i1] = wc{:};

% axis;
yc  = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
yic = skew_norm_pdf(x_values, u_ic, sigma_ic, lambda_ic);
yi1 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);
% y1 = alpha*yc + beta*yi1 + (1-alpha-beta)*yic;
y1  = w1c*yc + w1ic*yic + w1i1*yi1;
% y4 = skew_norm_pdf

% ax1 = subplot(2,1,1);
cla;
hold on;
plot(x_values,yc*w1c,'LineWidth',2);
plot(x_values,yic*w1ic,'LineWidth',2, 'LineStyle', '--');
plot(x_values,yi1*w1i1,'LineWidth',2);
plot(x_values,y1,'LineWidth',2);


histogram(s1,100,'Normalization','pdf', 'FaceColor', 'none');

legend({'dist\_correct'; 'dist\_half\_incorrect'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});


% legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'mixture2'; 'mixture3'; 'hist\_first'; 'hist\_second'});

hold off;

end


