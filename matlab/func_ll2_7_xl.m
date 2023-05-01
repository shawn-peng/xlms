
function ll = func_ll2_7_xl(S1, S2, ws, ...
    u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, ...
    u_c2, sigma_c2, lambda_c2, u_ic2, sigma_ic2, lambda_ic2, u_i12, sigma_i12, lambda_i12, ...
    u_i22, sigma_i22, lambda_i22)

%     lower_bound = 1e-27;
    
    P_c1 = skew_norm_pdf(S1, u_c, sigma_c, lambda_c);
    P_ic1 = skew_norm_pdf(S1, u_ic, sigma_ic, lambda_ic);
    P_i11 = skew_norm_pdf(S1, u_i1, sigma_i1, lambda_i1);
    P_c2 = skew_norm_pdf(S2, u_c2, sigma_c2, lambda_c2);
    P_ic2 = skew_norm_pdf(S2, u_ic2, sigma_ic2, lambda_ic2);
    P_i12 = skew_norm_pdf(S2, u_i12, sigma_i12, lambda_i12);
    P_i22 = skew_norm_pdf(S2, u_i22, sigma_i22, lambda_i22);
    
    pj = zeros(7, size(P_c1, 2));
%     p1 = (ws(1)+ws(2)) * P_c1 ...
%        + (ws(3)+ws(5)) * P_ic1 ...
%        + (ws(4)+ws(6)+ws(7)) * P_i11;
%     p2 = (ws(3)+ws(4)) * P_c2 ...
%        + (ws(1)+ws(6)) * P_ic2 ...
%        + (ws(2)+ws(5)) * P_i12 ...
%        + (ws(7)) * P_i22;
    pj(1, :) = ws(1) * P_c1  .* P_ic2;
    pj(2, :) = ws(2) * P_c1  .* P_i12;
    pj(3, :) = ws(3) * P_ic1 .* P_c2;
    pj(4, :) = ws(4) * P_i11 .* P_c2;
    pj(5, :) = ws(5) * P_ic1 .* P_i12;
    pj(6, :) = ws(6) * P_i11 .* P_ic2;
    pj(7, :) = ws(7) * P_i11 .* P_i22;
    
%     p1(p1 < lower_bound) = lower_bound;
%     p2(p2 < lower_bound) = lower_bound;
%     p1(p1==0) = min(p1(p1~=0));
%     p2(p2==0) = min(p2(p2~=0));
%     ll = mean(log(p1)) + mean(log(p2));
    ll = mean(log(sum(pj, 1)));
end
