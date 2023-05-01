function mat = generate_simulated_queryresult_xl_normal(M, ...
    alpha, beta, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_i1, sigma_i1, lambda_i1, ...
    u_i2, sigma_i2, lambda_i2)
%generate_simulated_queryresult Generate a simulated dataset
%
n = M;
% m = N;

c  = normrnd(u_c,  sigma_c, [n, 1]);
ic = normrnd(u_ic, sigma_ic, [n, 1]);
i1 = normrnd(u_i1, sigma_i1, [n, 1]);
i2 = normrnd(u_i2, sigma_i2, [n, 1]);

% [u, sigma, lambda] = sn_para_est(sort(c'));
% [u, sigma, lambda] = sn_para_est(sort(ic'));
% [u, sigma, lambda] = sn_para_est(sort(i1'));
% [u, sigma, lambda] = sn_para_est(sort(i2'));

figure;
hold on;
histogram(c, 'BinWidth', 1, 'Normalization', 'pdf')
histogram(ic, 'BinWidth', 1, 'Normalization', 'pdf')
histogram(i1, 'BinWidth', 1, 'Normalization', 'pdf')
histogram(i2, 'BinWidth', 1, 'Normalization', 'pdf')

z1 = rand(n, 1) < alpha;
z2 = rand(n, 1) < beta;

flag_c_ic = z1 .* z2;
flag_c_i1 = z1 .* (1-z2);

flag_ic_i1 = (1-z1) .* z2;
flag_i1_i2 = (1-z1) .* (1-z2);

su1 = flag_c_ic .* c  + flag_c_i1 .* c  + flag_ic_i1 .* ic + flag_i1_i2 .* i1;
su2 = flag_c_ic .* ic + flag_c_i1 .* i1 + flag_ic_i1 .* i1 + flag_i1_i2 .* i2;

SU = [su1'; su2'];

S = sort(SU, 'descend');

s1 = S(1,:);
s2 = S(2,:);

figure;
histogram(s1, 'BinWidth', 1, 'Normalization', 'pdf')
figure;
histogram(s2, 'BinWidth', 1, 'Normalization', 'pdf')


mat = [s1; s2];

end

