function y = quantizer_L_level(x, x_max, L);
delta = 2*x_max/L;
y = floor((x+x_max)/delta)*delta -x_max + delta/2;
for ind = 1: length(y)
    if y(ind) > (L/2-0.5)*delta;
        y(ind) = (L/2-0.5)*delta;
    end
    if y(ind) < -(L/2-0.5)*delta;
        y(ind) = -(L/2-0.5)*delta;
    end
end
