%% README
%This is attempt to directly construct equation system using the flatten
%and unflatten functions that relate grid positions to the "flattened"
%index in the temperature vector being solved for. 

%% CONFIG

nDefault = 40;
 %first two args give topmost righmost coord. Gives amount of vars
gradientBC = 1; %zero or one, depending on if one side should be insulated or not
nonLinearVersion = 0;  %zero or one, depending on if we solve nonlinear sys w/ fsolve
LineAlongBCPlot = 0;
tryFancyDelaunay = 0;
generateConvergencePlot = 1;
if generateConvergencePlot == 1
    nArray = 4:2:30;
    convergenceArr = zeros(1, length(nArray));
else
    nArray = nDefault;
end

convIdx = 1;

for n = nArray
%% System setup
h = 2/n;
maxIdx = flatten(n/2+1, n+1, n);
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
    BClogical = ones(length(A),1);
    for k = 1:length(A)
        [i,j] = unflatten(k,n);
        if isOnBoundary(i,j,n)
            BClogical(k) = 0;
        end
    end
    t = fsolve(@(t)(A*t+(t.*BClogical).^4-b'), ones(length(A),1));
end
    
temp = t.^2;
convergenceArr(convIdx) = sum(temp)*h^2;
convIdx = convIdx +1;
end

if generateConvergencePlot == 0
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
            x(k) = h*(i-1); y(k) = h*(j-1); z(k) = T(i,j);
            k = k+1;
        end 
    end 
end 

if tryFancyDelaunay == 1
%Cant get this to work... 
    figure
%     tri = delaunayTriangulation(x,y,z);
%     triplot(tri)

    %trisurf(DT.Points(IO,1),DT.Points(IO,2));
    Data = [x y z];
%     DT = delaunayTriangulation(x,y,z);
%     IO = isInterior(DT);
%     D = [x y z];
%     DT = delaunayTriangulation(D);
%     tri = delaunay(x,y);
    %tri = (isinterior(tri));

    %trisurf(tri,x,y,z)

else 
    tri = delaunay(x,y);
    trisurf(tri,x,y,z); 
end

%% Plot temperature between B and C
if LineAlongBCPlot == 1
    x = 1:n/2+1; x = sqrt(5) * x;
    y = zeros(1, n/2+1);
    for k = 1:n/2+1
        i = k;
        j = (k-1)*2+1;
        y(k) = T(i,j);
    end
    plot(x,y)
end 

else 
    plot(nArray, convergenceArr);
end
    






