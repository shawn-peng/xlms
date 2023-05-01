function [x1, x2, x3, x4] = unpack4(x)
    xc = num2cell(x);
    [x1, x2, x3, x4] = xc{:};
end