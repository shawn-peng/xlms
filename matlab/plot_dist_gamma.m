if exist('bootstrap_num')
    plot_dist_gamma_bootstrap
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

species_folder = [results_folder,species];
if ~ exist(species_folder)
    mkdir(species_folder)
    mkdir([species_folder,'/distplot'])
    mkdir([species_folder,'/fitcurv'])
    mkdir([species_folder,'/fdrcurv'])
end

figure('Position', [0,0,figw,figh]);
% figure;
hold on;
box on;

nsp = 1;
if (strcmp(method(1:5), '_3s4c'))
    nsp = 2;
end
% axis;
yc = normpdf(S1, u_c, sigma_c);
yi1 = gampdf(S1-gamma_i, a_i, 1/b_i);
y1 = alpha*yc + (1-alpha)*yi1;
% y4 = skew_norm_pdf
% subplot(nsp,2,1)
% cla;
% hold on;

l_hist_s1 = histogram(S1,'BinWidth',1,'Normalization','pdf','EdgeColor',[0.3010 0.7450 0.9330]);

l_pdf_c = plot(S1,yc*alpha,'Color',c_color,'LineWidth',linewidth);
l_pdf_i1 = plot(S1,yi1*(1-alpha),'Color',i1_color,'LineWidth',linewidth);
l_pdf_s1 = plot(S1,y1,'Color',mix_color,'LineWidth',linewidth);

xlim([0, xlim_high])

xlabel('-log(EValue)')
ylabel('PDF');
% legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'hist\_first'});

% saveas(gcf,[species_folder,'/distplot/',method,'.png'])

% figure('Position', [10,10,1920,1080]);
% figure;
% hold on;

% subplot(nsp, 2, 1);
% hold on;

m = size(S1,1);
h1emp = (1:m)' / m;

hc = normcdf(S1, u_c, sigma_c);
hi1 = gamcdf(S1-gamma_i, a_i, 1/b_i);
h1 = alpha*hc + (1-alpha)*hi1;

if plotcdf
    yyaxis('right');
    
    l_cdf_emp = plot(S1,h1emp, 'Color', cdfcolor);
    l_cdf_est = plot(S1,h1, 'Color', cdfcolor,'LineStyle',':');

    seq = flip(unique(omat(:,1)));

    n = size(seq,1);
    for i = 1:n
        s = seq(i);
        fdr = fdr_a(alpha, u_c, sigma_c, a_i, b_i, gamma_i, s);
    %     fdr
        if fdr > 0.1
            break
        end
    end
    % xline(s);
    % text(s,0.4,'\leftarrow 1%fdr')
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
% subplot(nsp, 2, 2);
% hold on;
% 
% S2 = sort(omat(:,2));
% m = size(S2,1);
% h2emp = (1:m) / m;
% pemp = plot(S2,h2emp);
% 
% hi2 = skew_norm_cdf(x_values, u_i2, sigma_i2, lambda_i2);
% h2 = alpha*hi1 + beta*hc + (1-alpha-beta)*hi2;
% pest = plot(x_values, h2);
% 
% legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});
% ylim([0,1]);

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


% saveas(gcf,[species_folder,'/fitcurv/',method,'.png'])
print([species_folder,'/distplot/',method,'.png'], '-dpng', '-r320'); 
print([species_folder,'/distplot/',method,'.eps'], '-depsc', '-r320'); 

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