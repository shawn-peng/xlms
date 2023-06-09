
function f = make_param_check_function_s2combined(x, ...
    w_11, w_12, w_21, w_22, ...
    u_1, sigma_1, lambda_1, ...
    u_2, sigma_2, lambda_2, ...
    constraints)

function debug()
    figure;
    hold on;
    plot_skewnorm(gca,x,u_1,sigma_1,lambda_1,w_11);
    plot_skewnorm(gca,x,u_2,sigma_2,lambda_2,w_12);
    yyaxis("right");
    plot_skewnorm_cdf(gca,x,u_1,sigma_1,lambda_1,w_11);
    plot_skewnorm_cdf(gca,x,u_2,sigma_2,lambda_2,w_12);
    figure;
    hold on;
    plot_skewnorm(gca,x,u_1,sigma_1,lambda_1,w_21);
    plot_skewnorm(gca,x,u_2,sigma_2,lambda_2,w_22);
end

check_pdf = false;
check_cdf = false;

for c = constraints
    if c == "pdf"
        check_pdf = true;
    elseif c == "cdf"
        check_cdf = true;
    end
end

function flag = helper(param_name, param_val, fuzzy)
    eval([param_name, '=param_val;']);
    
    mode_1 = skew_norm_mode(u_1, sigma_1, lambda_1);
    x_pdf = x(x>mode_1);
    x_cdf = x(:);


    if check_pdf
        flag = check_skewnorm_pdf_higher(x_pdf, ...
            w_11, w_12, ...
            u_1, sigma_1, lambda_1, ...
            u_2, sigma_2, lambda_2, ...
            fuzzy);
        if ~flag
            return
        end
    end
    if check_cdf
        flag = check_skewnorm_cdf_higher(x_cdf, ...
            u_1, sigma_1, lambda_1, ...
            u_2, sigma_2, lambda_2, ...
            fuzzy * 5e4);
        if ~flag
            return
        end
    end
    if check_pdf
        flag = check_skewnorm_pdf_higher(x_pdf, ...
            w_22, w_21, ...
            u_2, sigma_2, lambda_2, ...
            u_1, sigma_1, lambda_1, ...
            fuzzy);
        if ~flag
            return
        end
    end
    flag = true;

end

f = @helper;
end

