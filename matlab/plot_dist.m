if exist('bootstrap_num')
    plot_dist_bootstrap
    return
end

minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;


% s1 = omat(:,1);
% s1 = s1(s1~=0);
% s2 = omat(:,2);
% s2 = s2(s2~=0);
% s3 = omat(:,3);
% s3 = s3(s3~=0);

S1 = sort(omat(:,1));
S1 = S1(S1~=0);
S2 = sort(omat(:,2));
S2 = S2(S2~=0);
S3 = sort(omat(:,3));
S3 = S3(S3~=0);

% figw = 350;
% figh = 150;
figure('Position', [0,0,figw,figh]);
hold on;
box on;

% results_folder = 'test_search/fragger_results/';
% results_folder = 'test_search/est_results_nist/';
% results_folder = 'test_search/est_results/';

species_folder = [results_folder,species];

distplot_dir = [species_folder,'/distplot/'];
fitcurv_dir = [species_folder,'/fitcurv/'];
fdrcurv_dir = [species_folder,'/fdrcurv/'];
sddir = [species_folder,'/sdcdf/'];

if exist('bootstrap_num')
    distplot_dir = [distplot_dir,method,'/bootstrap/'];
    fitcurv_dir = [fitcurv_dir,method,'/bootstrap/'];
    fdrcurv_dir = [fdrcurv_dir,method,'/bootstrap/'];
    sddir = [sddir,method,'/bootstrap/'];
end

if ~exist(species_folder) || ~exist(distplot_dir)
    mkdir(species_folder)
    mkdir(distplot_dir)
    mkdir(fitcurv_dir)
    mkdir(fdrcurv_dir)
    mkdir(sddir)
end


nsp = 1;
if (strcmp(method(1:5), '_3s4c'))
    nsp = 2;
end
% axis;
yc = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
yi1 = skew_norm_pdf(S1, u_i, sigma_i, lambda_i);
y1 = alpha*yc + (1-alpha)*yi1;
% y4 = skew_norm_pdf
% subplot(nsp,2,1)
% cla;
% hold on;

l_hist_s1 = histogram(S1,'BinWidth',bin_width,'Normalization','pdf','EdgeColor',[0.3010 0.7450 0.9330]);

l_pdf_c = plot(S1,yc*alpha,'Color',c_color,'LineWidth',linewidth);
l_pdf_i1 = plot(S1,yi1*(1-alpha),'Color',i1_color,'LineWidth',linewidth);
l_pdf_s1 = plot(S1,y1,'Color',mix_color,'LineWidth',linewidth);

xlim([xlim_low, xlim_high])

xlabel('-log(EValue)')
ylabel('PDF');
% legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});


% if (strcmp(method(1:5), '_3s4c'))
%     subplot(nsp,2,3)
%     cla;
%     hold on;
%     yi3 = skew_norm_pdf(x_values, u_i3, sigma_i3, lambda_i3);
%     y3 = alpha*yi2 + (1-alpha)*yi3;
%     plot(x_values,alpha*yi2,'LineWidth',2);
%     plot(x_values,(1-alpha)*yi3,'LineWidth',2);
%     plot(x_values,y3,'LineWidth',2);
%     histogram(s3,100,'Normalization','pdf');
%     legend({'dist\_i2'; 'dist\_i3'; 'mixture'; 'hist\_third'});
% %     histogram(mat3(:,2),100,'Normalization','pdf');
% %     histogram(mat3(:,3),100,'Normalization','pdf');
% end

% saveas(gcf,[species_folder,'/distplot/',method,'.png'])

% figure('Position', [10,10,1920,1080]);
% figure;
% hold on;

% subplot(nsp, 2, 1);
% hold on;

% S1 = sort(omat(:,1));
m = size(S1,1);
h1emp = (1:m)' / m;

% hc = skew_norm_cdf(x_values, u_c, sigma_c, lambda_c);
% hi1 = skew_norm_cdf(x_values, u_i, sigma_i, lambda_i);
hc = skew_norm_cdf(S1, u_c, sigma_c, lambda_c);
hi1 = skew_norm_cdf(S1, u_i, sigma_i, lambda_i);
h1 = alpha*hc + (1-alpha)*hi1;

seq = flip(unique(omat(:,1)));
n = size(seq,1);

if plotcdf
    yyaxis('right');
    l_cdf_emp = plot(S1, h1emp, 'Color', cdfcolor);
    l_cdf_est = plot(S1, h1, 'Color', cdfcolor, 'LineStyle',':');

    for i = 1:n
        s = seq(i);
        fdr = fdr_x(alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
    %     fdr
        if fdr > 0.1
            break
        end
    end
    if plotthres
        xline(s);
        text(s,0.4,'\leftarrow 1%fdr')
    end
    sdcdf = calc_delta_cdf(S1, h1, h1emp);
    text(xlim_high-2, 0.8, ...
        sprintf('\\delta_{CDF} = %.2f', sdcdf), ...
        'HorizontalAlignment', 'right')
    xlabel('-log(EValue)');
    ylabel('CDF');
    ylim([0,1]);
    ax = gca;
    ax.YAxis(2).Color = cdfcolor;

    if legending
        legending = false;
        legend([l_hist_s1, l_pdf_c, l_pdf_i1, l_pdf_s1, l_cdf_emp, l_cdf_est],...
            {'top score'; 'correct'; 'incorrect'; 'mixture'; 'empirical CDF'; 'estimated CDF'});
    end
else
    if legending
        legending = false;
        legend([l_hist_s1, l_pdf_c, l_pdf_i1, l_pdf_s1],...
            {'top score'; 'correct'; 'incorrect'; 'mixture'});
    end
end
% legend([l_cdf_emp, l_cdf_est], {'empirical CDF'; 'estimated CDF'});

% subplot(nsp, 2, 2);
% hold on;

% if (strcmp(method(1:5), '_3s4c'))
%     subplot(nsp,2,3)
%     hold on;
%     
%     S3 = sort(omat(:,3));
%     m = size(S3,1);
%     h3emp = (1:m) / m;
%     pemp = plot(S3,h3emp);
% 
%     hi3 = skew_norm_cdf(x_values, u_i3, sigma_i3, lambda_i3);
%     h3 = (alpha+beta)*hi2 + (1-alpha-beta)*hi3;
%     pest = plot(x_values, h3);
% 
%     legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});
%     ylim([0,1]);
% end

% saveas(gcf,[species_folder,'/distplot/',method,'.png'])
print([species_folder,'/distplot/',method,'.png'], '-dpng', '-r320'); 
print([species_folder,'/distplot/',method,'.eps'], '-depsc', '-r320'); 

% saveas(gcf,[species_folder,'/fitcurv/',method,'.png'])
if (strcmp(method, '_2s3ci'))
    figure('Position', [0,0,figw,figh]);
    hold on;
    box on;
    
    yc2 = skew_norm_pdf(S2, u_c, sigma_c, lambda_c);
    yi21 = skew_norm_pdf(S2, u_i, sigma_i, lambda_i);
    yi22 = skew_norm_pdf(S2, u_i2, sigma_i2, lambda_i2);
    y2 = alpha*yi21 + beta*yc2 + (1-alpha-beta)*yi22;
    
    l_hist_s2 = histogram(S2,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.6,...
        'BinWidth',bin_width,'Normalization','pdf','EdgeColor',[0.3010 0.7450 0.9330]);
    
    l_pdf_c2 = plot(S2,yc2*beta,'Color',c_color,'LineWidth',linewidth);
    l_pdf_i21 = plot(S2,yi21*alpha,'Color',i1_color,'LineWidth',linewidth);
    l_pdf_i22 = plot(S2,yi22*(1-alpha-beta),'Color',i2_color,'LineWidth',linewidth);
    l_pdf_s2 = plot(S2,y2,'Color',mix_color,'LineWidth',linewidth);

    xlim([xlim_low, xlim_high])

    xlabel('-log(EValue)');
    ylabel('PDF');

    m = size(S2,1);
    h2emp = (1:m)' / m;

    hc2 = skew_norm_cdf(S2, u_c, sigma_c, lambda_c);
    hi21 = skew_norm_cdf(S2, u_i, sigma_i, lambda_i);
    hi22 = skew_norm_cdf(S2, u_i2, sigma_i2, lambda_i2);
    h2 = alpha*hi21 + beta*hc2 + (1-alpha-beta)*hi22;
    
    if plotcdf
        yyaxis('right');
        l_cdf_emp2 = plot(S2, h2emp, 'Color', cdfcolor);
        l_cdf_est2 = plot(S2, h2, 'Color', cdfcolor, 'LineStyle', ':');

        sdcdf = calc_delta_cdf(S2, h2, h2emp);
        text(xlim_high-2, 0.8, ...
            sprintf('\\delta_{CDF} = %.2f', sdcdf), ...
            'HorizontalAlignment', 'right')
        xlabel('-log(EValue)');
        ylabel('CDF');
        ylim([0,1]);
        ax = gca;
        ax.YAxis(2).Color = cdfcolor;
        
        if legending2
            legending2 = false;
            legend([l_hist_s2, l_pdf_c2, l_pdf_i21, l_pdf_i22, l_pdf_s2, l_cdf_emp2, l_cdf_est2],...
                {'second score'; 'correct'; 'incorrect 1'; 'incorrect 2'; 'mixture'; 'empirical CDF'; 'estimated CDF'});
        end
    else
        if legending2
            legending2 = false;
            legend([l_hist_s2, l_pdf_c2, l_pdf_i21, l_pdf_i22, l_pdf_s2],...
                {'second score'; 'correct'; 'incorrect 1'; 'incorrect 2'; 'mixture';});
        end
    end
    
%     saveas(gcf,[species_folder,'/distplot/',method,'_2.png'])
    print([species_folder,'/distplot/',method,'_2.png'], '-dpng', '-r320');
    print([species_folder,'/distplot/',method,'_2.eps'], '-depsc', '-r320'); 


end

sdcdf = calc_delta_cdf(S1, h1, h1emp);

sddir = [species_folder, '/sdcdf/'];
if ~exist(sddir)
    mkdir(sddir)
end
save([sddir, method, '.mat'], 'sdcdf');


function sdcdf = calc_delta_cdf(S1, h1, h1emp)
    dh = h1emp - h1;
    % sdh = diff(S1).*(1/2 x.^2)

    ddh = diff(dh);
    ds = diff(S1);
    k = ddh ./ ds;
    % k=(h1-h0)/(x1-x0);
    % s=integral[(h0+x*k)*dx]
    % =h0*x + 1/2 x.^2*k | x0=0, x1=ds
    n = size(dh,1);
    s = dh(1:n-1) .* ds + 1/2 * ds .* ddh;
    sdcdf = sum(abs(s));
end

% hc = skew_norm_cdf(S1, u_c, sigma_c, lambda_c);
% hi1 = skew_norm_cdf(S1, u_i, sigma_i, lambda_i);
% h1 = alpha*hc + (1-alpha)*hi1;


