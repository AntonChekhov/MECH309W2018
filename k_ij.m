function [val, partialx, partialy] = k_ij(x,y)
val = 2+cos(x+y);
partialx = -sin(x+y);
partialy = -sin(x+y);
end

