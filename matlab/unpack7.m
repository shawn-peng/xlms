function [x1, x2, x3, x4, x5, x6, x7] = unpack7(x)
    xc = num2cell(x);
    [x1, x2, x3, x4, x5, x6, x7] = xc{:};
end