
def param_binary_search(old_val, new_val, test_func, tolerance=1e-4):
    if test_func(new_val):
        return new_val
    if not test_func(old_val, 2e-7):
        print('old_val does not satisfy constraint, skip')
        return new_val
    assert test_func(old_val, 2e-7)

    while True:
        cur_val = (old_val + new_val) / 2
        if test_func(cur_val):
            old_val = cur_val
        else:
            new_val = cur_val
            cur_val = old_val

        if new_val == old_val or (new_val - old_val) / old_val < tolerance:
            break

    assert test_func(cur_val, 2e-7)
    return cur_val
