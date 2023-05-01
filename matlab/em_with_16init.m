
best_model = [];
best_init = [];
best_ll = -inf;
ll_list = {};
i = 1;
for sl0 = [1,-1]
    for sl1 = [1,-1]
        for sl2 = [1,-1]
            for sl3 = [1,-1]
                [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3] = EM3_3(mat3',sl0,sl1,sl2,sl3)
                S1 = mat3(:,1);
                S2 = mat3(:,2);
                S3 = mat3(:,3);
                ll = func_ll3_3(S1, S2, S3, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
                ll_list{i} = {ll, [sl0, sl1, sl2, sl3], [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3]};
                if ll > best_ll
                    best_ll = ll;
                    best_init = [sl0, sl1, sl2, sl3];
                    best_model = [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3];
                end
                i = i + 1;
            end
        end
    end
end

% [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = best_model
[alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3] = feval(@(x) x{:}, num2cell(best_model))

ll_list
