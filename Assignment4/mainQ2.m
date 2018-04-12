%% Setup
startx = 1; endx = 10;
N = 5; xpt = linspace(startx,endx,N+1); ypt = arrayfun(@(x) log(x), xpt); L = endx-startx; h = L/N;
A = zeros(3*N, 3*N);

rhs = [ypt zeros(1, 2*N-1)]'; 


A = [eye(N) zeros(N, 2*N);
    zeros(1, N-1) 1 zeros(1, N-1) h zeros(1, N-1) h^2 %interpol
    diag(ones(1,N-1)) + diag(-1*ones(1,N-2), 1) [zeros(N-2, 1);-1] h*eye(N-1) zeros(N-1, 1) h^2*eye(N-1) zeros(N-1, 1) %contrin
    zeros(N-1, N) diag(ones(1,N-1)) + diag(-1*ones(1,N-2), 1) [zeros(N-2, 1);-1] h*2*eye(N-1) zeros(N-1, 1)];

artificialConstraintA = [zeros(1, N) 1 zeros(1, 2*N-1)];
A = [A;artificialConstraintA];
%% Minimize length
min = 100000000; tmin = 0;
for t = linspace(-5, 5, 100)
    L=0; %I commit a sin by using fsolve but my gaussian elimination is being weird with this system and don't have time to debug
    rhs(length(rhs)) = t; sol = fsolve(@(x) A*x-rhs, ones(3*N, 1));
    arclength = 0;
    a = sol(1:N); b = sol(N+1:2*N); c = sol(2*N+1:3*N);
    for i = 1:N
        xtemp = linspace(xpt(i), xpt(i)+h, 10);
        Ltemp = arrayfun(@(x) sqrt(1+(b(i)+2*c(i)*x)^2), xtemp);
        Ltemp = sum(Ltemp)*h/10;
        L = L + Ltemp;
    end
    
    if L < min
        min = L; tmin = t;
    end
end
rhs(length(rhs)) = tmin; sol = gaussianElim(A, rhs); sol = fsolve(@(x) A*x-rhs, ones(3*N, 1));
%% Minimize area
areaMin = 100000000; tAmin = 0;
for t = linspace(-5, 5, 100)
    rhs(length(rhs)) = t; solA = fsolve(@(x) A*x-rhs, ones(3*N, 1));
    area = 0;
    a = solA(1:N); b = solA(N+1:2*N); c = solA(2*N+1:3*N);
    for i = 1:N
        xtemp = linspace(xpt(i), xpt(i)+h, 10);
        Atemp = arrayfun(@(x)  abs(a(i)+ b(i)*(x-xpt(i))+ c(i)*(x-xpt(i))^2), xtemp);
        Atemp = sum(Atemp)*h/10;
        area = area + Atemp;
    end
    
    if area < min
        areaMin = area; tAmin = t;
    end
end
rhs(length(rhs)) = tAmin; solA = gaussianElim(A, rhs); solA = fsolve(@(x) A*x-rhs, ones(3*N, 1));

x = []; y = []; yA =[];  a = sol(1:N); b = sol(N+1:2*N); c = sol(2*N+1:3*N);
%solutions for the area
aA = solA(1:N); bA = solA(N+1:2*N); cA = solA(2*N+1:3*N);

for i = 1:N
    xtemp = linspace(xpt(i), xpt(i)+h, 10);
    ytemp = arrayfun(@(x) a(i)+ b(i)*(x-xpt(i))+ c(i)*(x-xpt(i))^2, xtemp);
    yAtemp = arrayfun(@(x) aA(i)+ bA(i)*(x-xpt(i))+ cA(i)*(x-xpt(i))^2, xtemp);
    x = [x xtemp]; y = [y ytemp]; yA = [yA yAtemp]; 
end
plot(x,y); hold on
plot(x,yA); hold on
scatter(xpt, ypt); hold on
legend('Minimize length', 'Minimize area', 'Points considered')






























