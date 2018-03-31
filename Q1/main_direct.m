%% README
%This is attempt to directly construct equation system using the flatten
%and unflatten functions that relate grid positions to the "flattened"
%index in the temperature vector being solved for. 

%% CONFIG

n = 8;
h = 2/n;
maxIdx = flatten(n/2+1, n+1, n); %first two args give topmost righmost coord. Gives amount of vars
gradientBC = 0; %zero or one, depending on if one side should be insulated or not
nonLinearVersion = 1;  %zero or one, depending on if we solve nonlinear sys w/ fsolve


%% System setup
A = zeros(2, maxIdx);
t = zeros(1, maxIdx); b = zeros(1, maxIdx);
for k = 1:length(t)
    disp(k)
    [i,j] = unflatten(k,n);
    if ~isOnBoundary(i,j,n)
        %DISCRETIZED PDE CONSTRAINT
        flat_ip = flatten(i+1,j,n); flat_im = flatten(i-1,j,n); 
        flat_jp = flatten(i,j+1,n); flat_jm = flatten(i,j-1,n);
        flat_ij = k;
        [k_va, kpartialx, kpartialy] = k_ij((i-1)*h, (j-1)*h);
        f = f_ij((i-1)*h, (j-1)*h);
        
        A(k, flat_ij) = 4*k_va/h;
        A(k,flat_ip) = -(k_va/h+kpartialx);
        A(k,flat_im) = -(k_va/h-kpartialx);
        A(k,flat_jp) = -(k_va/h+kpartialy);
        A(k,flat_jm) = -(k_va/h-kpartialy);
        b(k) = 2*h*f;
    else 
        %SET BOUNDARY CONDITIONS
        if gradientBC == 0
            A(k,k) = 1; b(k) = 0;
        else 
            if isOnBoundary(i,j,n) == 1
                A(k,k) = 1; b(k) = 0;
            elseif isOnBoundary(i,j,n) == 2
                toto = flatten(i-1,j,n); 
                A(k,k) = 1; A(k,toto) = -1; 
                b(k) = 0;
            end
        end
    end 
   
end

%% Solve and plot
if nonLinearVersion == 0
    t = gaussianElim(A,b);
else
    t = fsolve(@(t)(NonLinearT4Func(A,t,b,n)), ones(length(A),1));
end

T = zeros(n+1, n+1);
for i = 1:n+1
    for j = n+1
        T(i,j) = missing;
    end
end
for k = 1:length(t)
    [i,j] = unflatten(k,n);
    T(i,j) = t(k);
end 

figure 
k = 1; 
x = zeros(length(t),1);
y = zeros(length(t),1);
z = t;
for i = 1:n+1 
    for j = 1:n+1
        if ~ismissing(T(i,j)) 
            x(k) = i; y(k) = j; z(k) = T(i,j);
            k = k+1;
        end 
    end 
end 

tri = delaunay(x,y);
%tri = (isinterior(tri));

trisurf(tri,x,y,z)





