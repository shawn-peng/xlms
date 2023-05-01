function plot_dist_gamma_fn(mat, alpha, u_c, sigma_c, a_i, b_i, gamma_i)
minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;

omat = mat;
s1 = omat(:,1);
s1 = s1(s1~=0);

% figure('Position', [10,10,1920,1080]);
cla;
hold on;

% axis;
yc = normpdf(x_values, u_c, sigma_c);
yi1 = gampdf(x_values-gamma_i, a_i, 1/b_i);
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


end
