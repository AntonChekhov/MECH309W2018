function [val, partialx, partialy] = k_ij(x,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
val = 2+cos(x+y);
partialx = -sin(x+y);
partialy = -sin(x+y);

% val = 1;
% partialx = 0;
% partialy = 0;

end

