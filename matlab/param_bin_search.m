function new_param = param_bin_search(old, new_try, cond_func)
tolerance = 1e-4;

if cond_func(new_try, 0)
    new_param = new_try;
    return;
end

assert(cond_func(old, 1e-7));
while true
    new_param = (old + new_try) / 2;
    if cond_func(new_param, 0)
        old = new_param;
    else
        new_try = new_param;
        new_param = old;
    end
    
    if new_try == old || (new_try - old) / old < tolerance
        break
    end
end
assert(cond_func(new_param, 1e-7));

end