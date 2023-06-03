dataset_name = 'Alban';

load(['../results/matdata/scoremats/', dataset_name, '.mat']);
mat2 = mat(:,mat(2,:)~=0);
best_ll = -inf;

s1 = mat2(1,:);
s1 = sort(s1, 'descend');
q1 = quantile(s1, 0.01);
s1 = s1(:, (s1 > q1));

info_file = sprintf('%s%s.json', '../results/info/', dataset_name);
info = load_json(info_file)

d = sprintf('../figures_diff_init/%s', dataset_name);
if ~exist(d)
    mkdir(d);
else
    delete([d, '/*']);
end

for c_range = 1/5:0.01:4/5
% for c_range = 0.62:0.01:4/5
    sl1 =  1;
    sl2 = -1;
    sl3 = -1;
    sl4 = -1;
    % [params, ll, ll1] = EM2_4ic_xl_diff_init(mat2, sl1, sl2, sl3, sl4, c_range);
    [params, ll, ll1] = EM2_4ic_xl( ...
        mat2, sl1, sl2, sl3, sl4, 'c_range', c_range);

    if ll > best_ll
        best_ll = ll;
        best_params = params;
        best_sls = [sl1, sl2, sl3, sl4];
    end
    
    subplot(2, 1, 1);
    
    ws = best_params{1};
    theta = best_params{2};
    [u_c, sigma_c, lambda_c] = unpack_skntheta(theta(1));
    [u_ic, sigma_ic, lambda_ic] = unpack_skntheta(theta(2));
    [u_i1, sigma_i1, lambda_i1] = unpack_skntheta(theta(3));
    [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta(4));
    
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
    plot(s1, fdr);
    
    ti = sprintf('%s\\_(%d,%d,%d,%d)_%.2f', escape_underscore(dataset_name), sl1, sl2, sl3, sl4, c_range);
    title(ti);
    filename = sprintf('%s/%d_%d_%d_%d_%.2f.png', d, sl1, sl2, sl3, sl4, c_range);
    saveas(gcf, filename);
end

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
ti = sprintf('%s\\_(%d,%d,%d,%d)', escape_underscore(dataset_name), sl1, sl2, sl3, sl4);

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
plot(s1, fdr);

%     plot(info.scores, info.curve_fdr)

title(ti);
filename = sprintf('%s/best.png', d);
saveas(gcf, filename);
close all;

function s = escape_underscore(s)
    s = strrep(s, '_', '\_');
end