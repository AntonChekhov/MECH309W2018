function [sum] = besselIntegral(y,n)
%function calculates the integral involved in the bessel function of degree n at a
%point y

%variables
intervSize = 10000;
dPhi = pi/intervSize; %gives subintervals
sum = 0;
for k = 1:intervSize  
    sum = sum +( (sin((dPhi/2) + (k-1)*dPhi))^(2*n) * cos( y*cos((dPhi/2) + (k-1)*dPhi) ) ); %evaluating the function and adding it to the total sum
end

sum = sum*dPhi; %multiply sum by differential to get the integral

end

