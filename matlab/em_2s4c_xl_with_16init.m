disp(['using 8 inits'])

best_model = [];
best_init = [];
best_ll = -inf;
ll_list = [];
for sl1 = [1,-1]
    for sl2 = [1,-1]
        for sl3 = [1,-1]
            for sl4 = [1,-1]
                sl1, sl2, sl3, sl4
                [alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2] = EM2_4_xl(mat',sl1,sl2,sl3,sl4);
                ll = func_ll2_1(mat', alpha, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2);
                ll_list = [ll_list, ll];
                if ll > best_ll
                    best_ll = ll;
                    best_init = [sl1, sl2];
                    best_model = [alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2];
                end
            end
        end
    end
end

[alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2] = feval(@(x) x{:}, num2cell(best_model))

ll_list
