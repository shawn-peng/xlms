method='_2s4ic_xl'

% species_folder = [results_folder,species]
dataset_res_folder = [results_dir, dataset_name]
dataset_fig_folder = [figures_dir, dataset_name]

% [alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = EM2_1(omat',1,1)
em_2s4ic_xl_with_16init

% plot_dist_xl
q1 = quantile(mat2(1,:), 0.01);
outlier_flags = mat2(1,:) <= q1;
mat2_new = mat2(:, ~outlier_flags);

figure('Position', [10,10,2000,500]);


ws = best_params{1};
theta = best_params{2};

[u_c, sigma_c, lambda_c] = unpack_skntheta(theta(1));
[u_ic, sigma_ic, lambda_ic] = unpack_skntheta(theta(2));
[u_i1, sigma_i1, lambda_i1] = unpack_skntheta(theta(3));
[u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta(4));

plot_dist_i_xl_fn(mat2, '_s2_c4_xl', ws, ...
    u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2);

sl1 = best_sls(1);
sl2 = best_sls(2);
sl3 = best_sls(3);
sl4 = best_sls(4);

if ~exist('bootstrap_num')
ti = sprintf('%s\\_(%d,%d,%d,%d)', ...
    escape_underscore(dataset_name), ...
    sl1, sl2, sl3, sl4);
else
ti = sprintf('%s\\_(%d,%d,%d,%d)\\_bs%03d', ...
    escape_underscore(dataset_name), ...
    sl1, sl2, sl3, sl4, ...
    bootstrap_num);
end

subplot(2, 1, 1);

fdr = fdr_xl(s1, ...
    ws(1), ws(2), ws(3), ...
    u_c, sigma_c, lambda_c, ...
    u_i1, sigma_i1, lambda_i1, ...
    u_ic, sigma_ic, lambda_ic);
sthres = fdr_thres(s1, fdr, 0.01);

xline(sthres)
text(sthres,0.003,'\leftarrow 1%fdr')

xline(info.fdr_thres)
text(info.fdr_thres,0.005,'\leftarrow TDA 1%fdr')

yyaxis right
hold on;
plot(s1, fdr);

plot(info.scores, info.curve_fdr, '--');

title(ti);
if ~exist('bootstrap_num')
figure_dir = sprintf('%s', figures_dir);
filename = sprintf('%s/%s.png', figure_dir, dataset_name);
else
figure_dir = sprintf('%s/%s/bootstrap', figures_dir, dataset_name);
filename = sprintf('%s/%03d.png', figure_dir, bootstrap_num);
end
if ~exist(figure_dir)
    mkdir(figure_dir)
end
saveas(gcf, filename);


if ~exist('bootstrap_num')
    param_folder = [dataset_res_folder,'/params/'];
    paramfile = [param_folder,method,'.mat'];
else
    param_folder = [dataset_res_folder,'/params/',method,'/bootstrap/'];
    paramfile = [param_folder,num2str(bootstrap_num),'.mat'];
end
if ~exist(param_folder)
    mkdir(param_folder)
end

save(paramfile, 'best_params');
% save(paramfile, ['alpha', 'u_c', 'sigma_c', 'lambda_c', 'u_i', 'sigma_i', 'lambda_i'])


