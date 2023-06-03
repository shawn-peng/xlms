function [x1, x2, x3, x4, x5, x6, x7, x8] = unpack8(x)
    xc = num2cell(x);
    [x1, x2, x3, x4, x5, x6, x7, x8] = xc{:};
end