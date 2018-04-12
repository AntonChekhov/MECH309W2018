function [val] = theta(x)
%function evaluates the theta function provided in the project

E=15000000000;
R=0.15;
S=pi*R^2;
I=(pi*R^4)/4;
rho=1000;
g=9.81;
q=rho*S*g/(E*I);
kappasqrd=4*q/9;

val = -1 * sqrt(x).* besselFn(sqrt(kappasqrd)*(x).^(3/2),-1/3);
end

