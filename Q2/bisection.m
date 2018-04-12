function [r] = bisection(a,b,eps,n)
%function performs bisection on the bessel function of order n.
m = 1;
a1 = a;
b1 = b;

while (b1 - a1) >= eps
    r = (a1+b1)/2;
    if besselFn(r,n)*besselFn(a1,n) > 0
        a1 = r;
    else
        b1 = r;
    end
    m = m + 1;
end

end

