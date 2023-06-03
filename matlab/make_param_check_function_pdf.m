
function f = make_param_check_function_pdf(x, ...
    w_1, w_2, ...
    u_1, sigma_1, lambda_1, ...
    u_2, sigma_2, lambda_2)

function debug()
    figure;
    hold on;
    plot_skewnorm(gca,x,u_1,sigma_1,lambda_1,w_1);
    plot_skewnorm(gca,x,u_2,sigma_2,lambda_2,w_2);
    yyaxis("right");
    plot_skewnorm_cdf(gca,x,u_1,sigma_1,lambda_1,w_1);
    plot_skewnorm_cdf(gca,x,u_2,sigma_2,lambda_2,w_2);
end

function flag = helper(param_name, param_val, fuzzy)
    eval([param_name, '=param_val;']);
    
    mode_1 = skew_norm_mode(u_1, sigma_1, lambda_1);
    x_pdf = x(x>mode_1);

    flag = check_skewnorm_pdf_higher(x_pdf, ...
        w_1, w_2, ...
        u_1, sigma_1, lambda_1, ...
        u_2, sigma_2, lambda_2, ...
        fuzzy);
end

f = @helper;
end

