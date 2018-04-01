function [val] = NonLinearT4Func(A,t,b,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

BClogical = ones(length(A),1);
for k = 1:length(A)
    [i,j] = unflatten(k,n);
    if isOnBoundary(i,j,n)
        BClogical(k) = 0;
    end
end

val = A*t+(t'.*BClogical).^4-b;
end