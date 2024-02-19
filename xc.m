function [x,F] = xc(t)
    x = cos(2*pi*1000*t.^3/3);
    F = 2*pi*(1000*t.^2);
end 