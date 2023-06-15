function [x1, x2, x3] = unpack3(x)
    xc = num2cell(x);
    [x1, x2, x3] = xc{:};
end