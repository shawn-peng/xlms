function mat = generate_simulated_queryresult_xl_nomax(M, ...
    alpha, beta, p_c_ic, p_c_i1, p_ic_i1, p_ic2_i1, ...
    u_c, sigma_c, lambda_c, ...
    u_ic, sigma_ic, lambda_ic, ...
    u_ic2, sigma_ic2, lambda_ic2, ...
    u_i1, sigma_i1, lambda_i1, ...
    u_i2, sigma_i2, lambda_i2)
%generate_simulated_queryresult Generate a simulated dataset
%
n = M;
% m = N;

w1 = alpha * beta * p_c_ic;
w2 = alpha * (1-beta) * p_c_i1;
w3 = alpha * beta * (1-p_c_ic);
w4 = alpha * (1-beta) * (1-p_c_i1);
w5 = (1-alpha) * beta * p_ic_i1 * (1-p_ic2_i1);
w6 = (1-alpha) * beta * (1-p_ic_i1);
w7 = (1-alpha) * beta * p_ic_i1 * p_ic2_i1;
w8 = (1-alpha) * (1-beta);

ws = [w1, w2, w3, w4, w5, w6, w7, w8];


c = randn_skew([n, 1], u_c, sigma_c, lambda_c);
ic = randn_skew([n, 1], u_ic, sigma_ic, lambda_ic);
i1 = randn_skew([n, 1], u_i1, sigma_i1, lambda_i1);
ic2 = randn_skew([n, 1], u_ic2, sigma_ic2, lambda_ic2);
i2 = randn_skew([n, 1], u_i2, sigma_i2, lambda_i2);

% [u, sigma, lambda] = sn_para_est(sort(c'));
% [u, sigma, lambda] = sn_para_est(sort(ic'));
% [u, sigma, lambda] = sn_para_est(sort(i1'));
% [u, sigma, lambda] = sn_para_est(sort(i2'));
% 
% figure;
% hold on;
% histogram(c, 'BinWidth', 1, 'Normalization', 'pdf')
% histogram(ic, 'BinWidth', 1, 'Normalization', 'pdf')
% histogram(i1, 'BinWidth', 1, 'Normalization', 'pdf')
% histogram(i2, 'BinWidth', 1, 'Normalization', 'pdf')

zr = rand(n, 1);
wsum = cumsum(ws);
wsum0 = wsum - ws;
z = zr < wsum & zr >= wsum0;

assert(sum(sum(z, 2) == 1));

z1 = z(:, 1);
z2 = z(:, 2);
z3 = z(:, 3);
z4 = z(:, 4);
z5 = z(:, 5);
z6 = z(:, 6);
z7 = z(:, 7);
z8 = z(:, 8);

su1 = (z1+z2) .* c + (z3+z5+z7) .* ic + (z4+z6+z8) .* i1;

su2 = (z3+z4) .* c + (z1+z6) .* ic + (z2+z5) .* i1 + z7 .* ic2 + z8 .* i2;

SU = [su1'; su2'];

S = SU;

s1 = S(1,:);
s2 = S(2,:);

figure;
histogram(s1, 'BinWidth', 1, 'Normalization', 'pdf')
figure;
histogram(s2, 'BinWidth', 1, 'Normalization', 'pdf')


mat = [s1; s2];

end

