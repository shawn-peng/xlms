
s1 = sort(omat(:,1));
m = size(s1,1);

% pdf

yc = alpha*skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
yi1 = (1-alpha)*skew_norm_pdf(s1, u_i, sigma_i, lambda_i);
y1 = yc + yi1;

% cdf

h1emp = (1:m)' / m;
% pemp = plot(s1, h1emp);

hc = alpha*skew_norm_cdf(s1, u_c, sigma_c, lambda_c);
hi1 = (1-alpha)*skew_norm_cdf(s1, u_i, sigma_i, lambda_i);
h1 = hc + hi1;
% pest = plot(s1, h1);

% fdr

n = size(s1,1);
flag = false;
flag01 = false;
flag10 = false;
fdr_curv = [];
for i = 1:n
    s = s1(i);
    fdr = fdr_x(alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
    if ~flag
        thres = s;
    end
    if fdr <= 0.01
        flag = true;
    end
    if ~flag01
        thres01 = s;
    end
    if fdr <= 0.001
        flag01 = true;
    end
    if ~flag10
        thres10 = s;
    end
    if fdr <= 0.001
        flag10 = true;
    end
    fdr_curv(i) = fdr;
%     fdr_curv(i,:) = [fdr, sum(s1>s), s];
end
fdr_curv = fdr_curv';

% deltacdf

dh = h1emp - h1;
% sdh = diff(S1).*(1/2 x.^2)

ddh = diff(dh);
ds = diff(s1);
% k = ddh ./ ds;
% k=(h1-h0)/(x1-x0);
% s=integral[(h0+x*k)*dx]
% =h0*x + 1/2 x.^2*k | x0=0, x1=ds
n = size(dh,1);
s = dh(1:n-1) .* ds + 1/2 * ds .* ddh;
sdcdf = sum(abs(s));


% ll
ll_1 = func_ll2_1(s1, alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i);


