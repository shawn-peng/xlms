function [m1, m2] = trunc_norm_moments(U, sigma)
%trunc_norm_pdf pdf of truncated norm distribution
%   Detailed explanation goes here
    cdf = normcdf(U/sigma);
    flags = cdf == 0;
    
%     cdf(flags) = min(min(cdf(~flags)));
    pdf = normpdf(U/sigma);
    p = pdf ./ cdf;
    p(flags) = abs(U(flags)/sigma);
    m1 = U + sigma * p;
    m2 = U.^2 + sigma^2 + sigma * U .* p;
%     m2(isnan(m2)) = inf;
end

