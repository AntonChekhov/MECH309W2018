function [val] = integrand(n,t)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
val = t.^(n - (1/2)) .* exp(-t);
end

