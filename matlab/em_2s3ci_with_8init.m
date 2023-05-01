
S = omat';
s1 = S(1,:);
s1 = s1(s1~=0);
s2 = S(2,:);
s2 = s2(s2~=0);

best_model = [];
best_init = [];
best_ll = -inf;
ll_list = [];
for sl1 = [1,-1]
    for sl2 = [1,-1]
        for sl3 = [1,-1]
            [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2] = EM3(omat',sl1,sl2,sl3)
            
            ll = func_ll3(s1, s2, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2)
            ll_list = [ll_list, ll];
            if ll > best_ll
                best_ll = ll;
                best_init = [sl1, sl2, sl3];
                best_model = [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2];
            end
        end
    end
end

% [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = best_model
[alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2] = feval(@(x) x{:}, num2cell(best_model))

ll_list
