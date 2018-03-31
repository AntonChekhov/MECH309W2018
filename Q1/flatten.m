function [k] = flatten(i,j,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
botRectangleNodes = (n+1)*(n/2+1);
if (i > n/2+1 && j>n/2+1)
    k = missing;
    disp("Indices are outside of plate. WRONG!")
elseif j<= n/2 + 1
    k = (j-1)*(n+1) + i;
else
    tempj = j - (n/2+1);
    tempk = (tempj-1)*(n/2+1) + i;
    k = botRectangleNodes + tempk; %(i-2 - n/2)*(n/2+1) + j;
end
end

