function [best_theta, best_ll] = param_rand_search(theta_old, theta_new, cond_func, ll_func, sample_size)
    theta_range = [theta_old; theta_new];
    theta_range = sort(theta_range);

    best_ll = -inf;
    best_theta = [];
    assert(cond_func(theta_old))
    if cond_func(theta_new)
        best_theta = theta_new;
        best_ll = ll_func(theta_new);
        return;
    end
    for i = 1:sample_size
        theta = rand_theta(theta_range);
        cond = cond_func(theta);
        if cond
            ll = ll_func(theta);
            if ll > best_ll
                best_ll = ll;
                best_theta = theta;
            end
        end
    end
    ntries = 0;
    while numel(best_theta) == 0
        theta = rand_theta(theta_range);
        cond = cond_func(theta);
        if cond
            ll = ll_func(theta);
            if ll > best_ll
                best_ll = ll;
                best_theta = theta;
            end
            break
        end
        ntries = ntries + 1;
        if ntries > 1e4
            theta
        end
    end
        
end

function theta = rand_theta(theta_range)
    m = size(theta_range, 2);
    theta = rand([1,m]) .* diff(theta_range) + theta_range(1,:);
end