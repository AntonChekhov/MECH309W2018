function [val] = gammaFn(y)
%calculates gamma(y) where y = n + 1/2
if (y>1)
    val = integralSum(0,300,100000,y);
else
    val = (1/y) * gamma(y+1);
end

