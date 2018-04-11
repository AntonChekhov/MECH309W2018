%% Setup
startx = 1; endx = 10;
N = 5; xpt = linspace(startx,endx,N+1); ypt = arrayfun(@(x) log(x), xpt); L = xpt(length(xpt)); h = L/N;
c = ones(3*N, 1);
A = zeros(3*N, 3*N);

rhs = [ypt zeros(1, 2*N-1)]'; 

Atop = [diag(ones(1,length(ypt)-1)) zeros(N, 2*N);
    zeros(1, N-1) 1 zeros(1, N-1) h zeros(1, N-1) h^2;
    diag(ones(1,N-1)) + diag(-1*ones(1,N-2), 1) [zeros(N-2, 1);-1] h*eye(N-1) zeros(N-1, 1) h^2*eye(N-1) zeros(N-1, 1);
    zeros(N-1, N) diag(ones(1,N-1)) + diag(-1*ones(1,N-2), 1) [zeros(N-2, 1);-1] h^2*eye(N-1) zeros(N-1, 1);
    ];

BClength = [zeros(1, N-1) 1, zeros(1, N-1) 1 zeros(1, N-1) 1];
BCArea = [1/2 zeros(1, N-1) h/3 zeros(1, N-1) h^2/4 zeros(1, N-1)];

A = [Atop;BClength];
A = [Atop;BCArea];

sol = gaussianElim(A, rhs);

a = sol(1:N); b = sol(N+1:2*N); c = sol(2*N+1:3*N);

x = []; y = [];
i = 1;
for i = 1:N
    
    xtemp = linspace(xpt(i), xpt(i)+h, 10);
    ytemp = a(i) + b(i)*(xtemp-xpt(i)) + c(i)*(xtemp-xpt(i)).^2;
    x = [x xtemp];
    y = [y ytemp];
 
    
end

%% test
% i = 1;
% f = @(xtemp) (a(i) + b(i)*xtemp + c(i)*xtemp.^2);
% plot(x, arrayfun(@(x) f(x), x))




%% Plot
figure
plot(x,y); hold on
scatter(xpt, ypt)
