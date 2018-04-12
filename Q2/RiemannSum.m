function [sum] = RiemannSum(a,b,N)
%The function gives the integral of theta(x) from a to b by dividing the domain
%into N subintervals

%variables
dX = (b-a)/N; %gives subintervals
sum = 0;

for i = 1:N  
    sum = sum + theta(a + (dX/2) + (i-1)*dX); %evaluating f and adding it to the total sum
end

sum = sum.*dX; %multiply sum by differential to get the integral

end