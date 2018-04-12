%% Eigenvalues
f = @(x) (1-x)*((4-x)*(6-x)-1)-2*(2*(6-x)+3)-3*(2+3*(4-x));
x = linspace(-10, 10, 1000);
y = arrayfun(@(x) f(x), x);
plot(x,y); hold on
plot(x, zeros(1,length(x)))
r1 = fsolve(@(x)                     f(x), -5);
r2 = fsolve(@(x) f(x), 10);
r3 = fsolve(@(x) f(x), 5);                 
%% Eigenvectors
A = [1 2 -3; 2 4 1;-3 1 6];   
%Am very fsolve happy this assignment cause my gaussian elimination keeps
%breaking and aint no one got time fo dat. 

v1 = fsolve(@(x) (A-r1*eye(3))*x, ones(3,1));
v2 = fsolve(@(x) (A-r2*eye(3))*x, ones(3,1));
v3 = fsolve(@(x) (A-r3*eye(3))*x, ones(3,1));

%% 
A = [0 1 1; 0 1 1;1 0 0];   
r1 = (1+sqrt(5))/2; r2 = (1-sqrt(5))/2;
v1 = fsolve(@(x) (A-r1*eye(3))*x, ones(3,1));
v2 = fsolve(@(x) (A-r2*eye(3))*x, ones(3,1));


                                                        