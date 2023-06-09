
datasets = {
    % 'Alinden';
    'ecoli_xl';
    'MS2000225';
    'alban';
    'RPA';
    'CPSF';
    'ALott';
    'KKT4';
    'QE';
    'D1810';
    'peplib';
    'Gordon';
};

fails = {};

for i = 1:size(datasets)
    dataset_name = datasets{i};
    load(['../results/matdata/scoremats/', dataset_name, '.mat']);
    mat2 = mat(:,mat(2,:)~=0);
    best_ll = -inf;

    s1 = mat2(1,:);
    s1 = sort(s1, 'descend');
    q1 = quantile(s1, 0.01);
    s1 = s1(:, (s1 > q1));

    info_file = sprintf('%s%s.json', '../results/info/', dataset_name);
    info = load_json(info_file)

    d = sprintf('../figures_less_c2_ic2/%s', dataset_name);
    if ~exist(d)
        mkdir(d);
    else
        delete([d, '/*']);
    end

%         [params, ll, ll1] = EM2_4ic_xl(mat2,0, 0,0,0);
    for j = 0:3^5-1
        sl1 = get_sign(j, 1);
        sl2 = get_sign(j, 2);
        sl3 = get_sign(j, 3);
        sl4 = get_sign(j, 4);
        sl5 = get_sign(j, 5);
        if sl1 < 0 || sl2 < 0 || sl3 ~= sl4
            continue
        end

        ti = sprintf('%s\\_(%d,%d,%d,%d,%d)', escape_underscore(dataset_name), sl1, sl2, sl3, sl4, sl5);
        
        try
            % [params, ll, ll1] = EM2_5ic_xl_less_c2_2ic_rands(mat2, ...
            %     sl1,sl2,sl3,sl4,sl5,'tolerance',1e-4,'c_range',0.3, ...
            %     'constraints', []);
            % [params, ll, ll1] = EM2_5ic_xl_less_c2_2ic(mat2, ...
            %     sl1,sl2,sl3,sl4,sl5,'tolerance',1e-4,'c_range',0.3, ...
            %     'constraints', ["pdf"]);
            [params, ll, ll1] = EM2_5ic_xl_less_c2_2ic_rands(mat2, ...
                sl1,sl2,sl3,sl4,sl5,'tolerance',1e-4,'c_range',0.3, ...
                'constraints', ["pdf" "cdf"]);
            % [params, ll, ll1] = EM2_5ic_xl_less_c2_2ic_rands(mat2, ...
            %     sl1,sl2,sl3,sl4,sl5,'tolerance',1e-4,'c_range',0.3, ...
            %     'constraints', ["cdf"]);
        catch ME
            fails{end+1} = ti;
            continue
        end
        if ll > best_ll
            best_ll = ll;
            best_params = params;
            best_sls = [sl1, sl2, sl3, sl4, sl5];
        end
        
        subplot(2, 1, 1);
        
        ws = params{1};
        theta = params{2};
        [u_c, sigma_c, lambda_c] = unpack_skntheta(theta(1));
        [u_ic, sigma_ic, lambda_ic] = unpack_skntheta(theta(2));
        [u_ic2, sigma_ic2, lambda_ic2] = unpack_skntheta(theta(3));
        [u_i1, sigma_i1, lambda_i1] = unpack_skntheta(theta(4));
        [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta(5));
        
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
        
        title(ti);
        filename = sprintf('%s/%d_%d_%d_%d_%d.png', d, sl1, sl2, sl3, sl4, sl5);
        saveas(gcf, filename);
    end
    
    figure('Position', [10,10,2000,500]);
    
    ws = best_params{1};
    theta = best_params{2};
    [u_c, sigma_c, lambda_c] = unpack_skntheta(theta(1));
    [u_ic, sigma_ic, lambda_ic] = unpack_skntheta(theta(2));
    [u_ic2, sigma_ic2, lambda_ic2] = unpack_skntheta(theta(3));
    [u_i1, sigma_i1, lambda_i1] = unpack_skntheta(theta(4));
    [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta(5));
    plot_dist_i_xl_ic2_fn(mat2, '_s2_c4_xl', ws, ...
        u_c, sigma_c, lambda_c, ...
        u_ic, sigma_ic, lambda_ic, ...
        u_ic2, sigma_ic2, lambda_ic2, ...
        u_i1, sigma_i1, lambda_i1, ...
        u_i2, sigma_i2, lambda_i2);
    
    sl1 = best_sls(1);
    sl2 = best_sls(2);
    sl3 = best_sls(3);
    sl4 = best_sls(4);
    sl5 = best_sls(5);
    ti = sprintf('%s\\_(%d,%d,%d,%d,%d)', escape_underscore(dataset_name), sl1, sl2, sl3, sl4, sl5);
    
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
%     EM1_3i_xl(mat, 1, -1, -1)
%     EM2_4_xl(mat,1,-1,-1,-1)
%     EM2_7j_xl(mat2,1,1,-1,-1)
%     dlls = EM2_7uj_xl_new(mat2,1,1,-1,-1);
%     dlls = EM2_7j_xl_new(mat2,1,1,-1,-1);
%     dlls = EM2_7j_xl_new2(mat2,1,1,-1,-1);
%     dlls = EM2_10j_xl(mat2,[1,-1,-1,-1]);
%     dlls = EM2_4i_xl(mat2,1,-1,-1,-1);
%     dlls = EM2_4ic_xl(mat2,1,-1,-1,-1);
%     dlls = EM2_4ic_xl(mat2,1,1,1,-1);
%     EM3(mat,1,-1,-1)
%     figure;
%     plot(dlls)
end

function s = get_sign(j, d)
d = d - 1;
j = floorDiv(j, 3^d);
j = mod(j, 3);
if j == 0
    s = 0;
elseif j == 1
    s = -1;
else
    s = 1;
end
end

function s = escape_underscore(s)
    s = strrep(s, '_', '\_');
end
