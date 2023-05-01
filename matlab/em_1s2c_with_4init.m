disp(['using 4 inits'])

best_model = [];
best_init = [];
best_ll = -inf;
ll_list = [];
for sl1 = [1,-1]
    for sl2 = [1,-1]
		sl1, sl2
        [alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = EM2_1(mat',sl1,sl2);
        ll = func_ll2_1(mat', alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i);
        ll_list = [ll_list, ll];
        if ll > best_ll
            best_ll = ll;
            best_init = [sl1, sl2];
            best_model = [alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i];
        end
    end
end

% [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = best_model
[alpha, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = feval(@(x) x{:}, num2cell(best_model))

ll_list
