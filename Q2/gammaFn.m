function [val] = gammaFn(y)
%calculates gamma(y) where y = n + 1/2, using the Lanczos approximation

p0=1.000000000190015;
p1=76.18009172947146;
p2=-86.50532032941677;
p3=24.01409824083091;
p4=-1.231739572450155;
p5=1.208650973866179 * 10^(-3);
p6=-5.395239384953 * 10^(-6);

firstTerm = (sqrt(2*pi)) / y;
secTerm =  p0 + (p1/(y+1))+ (p2/(y+2))+ (p3/(y+3))+ (p4/(y+4))+ (p5/(y+5))+ (p6/(y+6));
thirTerm = (y+5.5)^(y+0.5);
fourTerm = exp(-y-5.5);
val = firstTerm*secTerm*thirTerm*fourTerm;
end
