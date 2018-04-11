y_nod = -1.5;
h = 0.25;
yEuler = ones(1, 4/h); yT = ones(1, 4/h);
yEuler(1) = y_nod; yT(1) = y_nod;
k = 2;
for t = h:h:4
    yEuler(k) = yEuler(k-1)+h*yprime(t, yEuler(k-1));
    yT(k)=yT(k-1)+h/2*(yprime(t,yT(k-1))+yprime(t+h,yT(k-1)+h*yprime(t,yT(k-1))));

    k = k+1;
end
t = 0:h:4;

figure
plot(t, yEuler)
hold on
plot(t, yT);
%% First piece analytic


f = @(t)(1/2*(cos(t)+sin(t))- 2*exp(-t));
plot(t, arrayfun(@(t) f(t), t));
hold on

firstzero=fsolve(@(t) f(t), 1);

%% Second piece
T = 1.0854;
ICsolver = @(a)( ...
                -cos(T)*0.5+sin(T)*0.5+ ... 
                a*exp(T));
A = fsolve(@(a) ICsolver(a), 1);

f2 = @(t)( ...
                -cos(t)*0.5+sin(t)*0.5+ ... 
                A*exp(t)); f2(T)
plot(t, arrayfun(@(t) f2(t), t));
secondzero = fsolve(@(t) f2(t), 5);
hold on

%% Third piece 
T = secondzero;
ICsolver = @(a)( ...
                cos(T)*0.5+sin(T)*0.5+ ... 
                a*exp(-T));
A = fsolve(@(a) ICsolver(a), 1);
f3 = @(t)( ...
                cos(t)*0.5+sin(t)*0.5+ ... 
                A*exp(-t)); f3(T)
plot(t, arrayfun(@(t) f3(t), t)); hold on
thirdzero = fsolve(@(t) f(t), 1);
hold on 
legend('f1', 'f2', 'f3')
plot(t, zeros(length(t),1))

%% Whole plot
%Y is analytical solution
figure
Y = ones(1,length(t));
Y(t <= firstzero) = arrayfun(@(t) f(t), t(t <= firstzero));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   arrayfun(@(t) f(t),  t(t <= firstzero));
Y(t >= firstzero & t <= secondzero) = arrayfun(@(t) f2(t),  t(t >= firstzero & t <= secondzero));
Y(t >= secondzero & t <= 4) = arrayfun(@(t) f3(t),  t(t >= secondzero & t <= 4));

plot(t, Y)
hold on
plot(t, arrayfun(@(t) f(t), t));
legend('Piecewise f1, f2, f3', 'f1 all the way')

%% Numerical vs analytical 
figure
plot(t,yEuler); hold on
plot(t,yT); hold on
plot(t, Y, '-o') 

legend('Euler method','Explicit Trapezoid', 'Quasi-exact')
ylabel('y')
xlabel('t')

%% Euler convergence
y_nod = -1.5;
nspace = 20:1000;
error = zeros(1,length(nspace));
figure
for i = 1:length(nspace)
    n = nspace(i); h = 1/n;
    yEuler = ones(1, n); yEuler(1) = y_nod; 
    k = 2;
    for t = h:h:4
        yEuler(k) = yEuler(k-1)+h*yprime(t, yEuler(k-1));
        k = k+1;
    end
    
    t = linspace(0,4,n);
    Y = ones(1,length(t));
    Y(t <= firstzero) = arrayfun(@(t) f(t), t(t <= firstzero));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   arrayfun(@(t) f(t),  t(t <= firstzero));
    Y(t >= firstzero & t <= secondzero) = arrayfun(@(t) f2(t),  t(t >= firstzero & t <= secondzero));
    Y(t >= secondzero & t <= 4) = arrayfun(@(t) f3(t),  t(t >= secondzero & t <= 4));
    
    error(i) = abs( ...
                sqrt(sum(Y.^2)) - sqrt(sum(yEuler.^2)) ...
                )/n;
    
end 


%% Log slope compute
figure
hspace = 4 ./ nspace;

loglog(hspace , error) 
errp = abs(error(2:length(error)));
errm = abs(error(1:(length(error)-1)));

plot(errp, errm)
% 
% slope = ( ...
%         (log(errdiff(900))-log(errdiff(700)))/ ...
%         (log(hspace(900))-log(hspace(700)))...
%         );





