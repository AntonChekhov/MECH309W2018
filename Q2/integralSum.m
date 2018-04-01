function [sum] = integralSum(a,b,intervSize,n)
%The function gives the integral of f from a to b by dividing the domain
%into N subintervals

%variables
dX = (b-a)/intervSize; %gives subintervals
sum = 0;

for i = 1:intervSize  
    sum = sum + integrand(n,a + (dX/2) + (i-1)*dX); %evaluating f and adding it to the total sum
end

sum = sum*dX; %multiply sum by differential to get the integral

end