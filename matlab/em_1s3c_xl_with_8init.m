disp(['using 8 inits'])

best_model = [];
best_init = [];
best_ll = -inf;
ll_list = [];
for sl1 = [1,-1]
    for sl2 = [1,-1]
        for sl3 = [1,-1]
            sl1, sl2, sl3
            [alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i, sigma_i, lambda_i] = EM1_3_xl(mat',sl1,sl2,sl3);
            ll = func_ll1_3_xl(mat(:,1)', alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_ic, sigma_ic, lambda_ic);
            ll_list = [ll_list, ll];
            if ll > best_ll
                best_ll = ll;
                best_init = [sl1, sl2];
                best_model = [alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i, sigma_i, lambda_i];
            end
        end
    end
end

[alpha, beta, u_c, sigma_c, lambda_c, u_ic, sigma_ic, lambda_ic, u_i, sigma_i, lambda_i] = feval(@(x) x{:}, num2cell(best_model))

ll_list
