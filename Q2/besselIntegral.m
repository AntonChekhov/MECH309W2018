function [sum] = besselIntegral(y,n)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
intervSize = 10000;
%variables
dPhi = pi/intervSize; %gives subintervals
sum = 0;
for k = 1:intervSize  
    sum = sum + cos(y*cos((dPhi/2) + (k-1)*dPhi))*(sin((dPhi/2) + (k-1)*dPhi))^(2*n); %evaluating f and adding it to the total sum
end

sum = sum*dPhi; %multiply sum by differential to get the integral

end

