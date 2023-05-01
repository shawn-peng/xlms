syms sym_alpha sym_beta sym_gamma sym_delta sym_epsilon eta1 eta2 eta3
syms C1 I11 C2 I21 I22 I32 I33
syms sym_M1 sym_M2 sym_M3

% Ab = [
%     %alpha, beta,   gamma,  delta,  epsi,   et1,    eta2,   eta3,   b
%     C1+I21, 0,      0,      0,      0,      1,      0,      -1,     0;
%     0,      C2,     0,      0,      0,      1,      -1,     -1,     0;
%     0,      0,      I22+I33,0,      0,      1,      -1,     0,      0;
%     0,      0,      0,      I11,    0,      0,      1,      0,      0;
%     0,      0,      0,      0,      I32,    0,      0,      1,      0;
%     1,      1,      1,      0,      0,      0,      0,      0,      1;
%     1,      1,      1,      0,      0,      0,      0,      0,      1;
%     1,      1,      1,      0,      0,      0,      0,      0,      1;
%     ]

I11 = sym_M1 - C1;
I22 = sym_M2 - C2 - I21;
I33 = sym_M3 - I32;
% assume(C1 + I11 == M1);
% assume(C2 + I21 + I22 == M2);
% assume(I32 + I33 == M3);
% assume(sym_alpha > 0);
% assume(sym_beta > 0);
% assume(sym_gamma > 0);
% assume(sym_delta > 0);
% assume(sym_epsilon > 0);

eq1 = 1/sym_alpha * (C1 + I21) + eta1 - eta3 == 0;
eq2 = 1/sym_beta * C2 + eta1 - eta2 - eta3 == 0;
eq3 = 1/sym_gamma * (I22 + I33) + eta1 - eta2 == 0;
eq4 = 1/sym_delta * I11 + eta2 == 0;
eq5 = 1/sym_epsilon * I32 + eta3 == 0;
eq6 = sym_alpha + sym_beta + sym_gamma == 1;
eq7 = sym_beta + sym_gamma == sym_delta;
eq8 = sym_alpha + sym_beta == sym_epsilon;

weights_eqs = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8, ...
    sym_alpha, sym_beta, sym_gamma, sym_delta, sym_epsilon, eta1, eta2, eta3);
