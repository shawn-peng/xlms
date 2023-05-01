function [u,sigma,lambda] = sn_para_est(X)
%sn_para_est estimate the parameters of the skew normal with method of
%moments
%   Detailed explanation goes here

m1 = mean(X);
m2 = var(X);
m3 = skewness(X);

syms sym_gamma3 sym_delta;

eqn = m3 == (4-pi)*2 * (sym_delta * sqrt(2/pi)) / (1-2*sym_delta^2/pi)^(2/3);
delta = vpasolve(eqn, sym_delta);
delta = double(delta);

lambda = delta / sqrt(1 - delta^2);

a1 = sqrt(2/pi);
b1 = (4/pi - 1) / a1;

% delta = sign(m3) ./ sqrt(a1^2 + m2 .* (b1 / abs(m3)).^(2/3));

sigma = sqrt(m2 ./ (1 - a1^2 * delta.^2));

u = m1 - a1*delta*sigma;

u_c = u;
sigma_c = sigma;
lambda_c = lambda;

n = size(X, 2);

s1 = X;

prev_ll = 0;
ll = -1;

figure;
hold on;
while abs(ll - prev_ll) > 1e-5
    prev_ll = ll;

    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);

    [Vc1 , Wc1 ] = trunc_norm_moments(delta_c  / sigma_c  * (s1-u_c ), sqrt(1-delta_c ^2));

    u_c_new  = sum( (s1 - Vc1  * Delta_c ) ) / n;

    Delta_c_new  = sum( Vc1  .* (s1 - u_c_new) ) / sum(Wc1);

    Gamma_c_new  = sum( (s1 - u_c_new).^2 - ...
                         2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + ...
                         Wc1 * Delta_c_new^2 ) / n ;

    u_c = u_c_new;
    Delta_c = Delta_c_new;
    Gamma_c = Gamma_c_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    
    p1 = skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
    cla;
    plot(s1, p1);
    histogram(s1, 'BinWidth', 1, 'Normalization', 'pdf');
    pause(.001);
    ll = mean(log(p1))
%     break
end

close


u = u_c;
sigma = sigma_c;
lambda = lambda_c;

end

