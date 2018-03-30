function [i,j] = unflatten(k, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

botRectangleNodes = (n+1)*(n/2+1);
j = floor((k-1)/(n+1))+1;

if j <= n/2 + 1
    i = k - (n+1)*(j-1);
else 
    tempK = k - botRectangleNodes;
    j = floor((tempK-1)/(n/2+1))+1;
    i = tempK - (n/2+1)*(j-1);
    j = j + n/2 + 1;
end 

if (i > n+1 || j>n+1)
    i = missing; j = missing;
    disp("Indices are outside of plate. WRONG!")
end
    




end

