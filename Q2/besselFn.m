function [val] = besselFn(y,n)

if n>=0
    val =  ((y).^n) /( (2^n) * (sqrt(pi))  * (gammaFn(n+1/2)));
    val = val .* besselIntegral(y,n);
else
    val = (2*(n+1)./y) .* besselFn(y,n+1) - besselFn(y,n+2); %recursive relation in case n is negative.
end


end

