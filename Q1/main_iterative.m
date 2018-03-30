%%  README

%Instead of setting up and solving equation system, use iterative approach
%Iterative approach idea given by Guangwei Liu.
%Algebraic proof of convergence is still required. 

n = 16;
L = 2;
tol = 10^-4;
gradientBC = 0; %zero or one, depending on if one side should be insulated or not
nonLinearVersion = 0;  %zero or one, depending on if we solve nonlinear sys w/ fsolve
%Define matrices for parameters that need to be evaluated everywhere

h = L/n;
T = ones(n+1, n+1); temp = zeros(n+1, n+1); residual = 1; maxIter = 10000; 

for iter = 1:maxIter
    if residual > tol
        temp = T;
        
        for i = 1:n+1
            for j = 1:n+1
                [k, kpartialx, kpartialy] = k_ij((i-1)*h, (j-1)*h);
                f = f_ij((i-1)*h, (j-1)*h);
                    if gradientBC == 1
                        if i == 1 || j == 1 || ... %Left and bottom sides
                             (i == n+1 && j <= n/2) || ... %Right and top sides
                             (i <= n/2 && j == n+1) || ... 
                             (i >= n/2 && j == n/2 ) %Inward corner sides
                             T(i,j) = 0;
                        elseif (i == n/2 && j >= n/2 )
                             T(i,j) = T(i-1, j);
                        elseif i > n/2 && j > n/2 
                             T(i,j) = missing;
                        else 
                            if nonLinearVersion == 0
                                 T(i,j) = h/(4*k)*...
                                    ( ...
                                        -2*h*f + ...
                                        ( ...
                                            T(i+1,j) * (k/h+kpartialx) + T(i-1,j) * (k/h-kpartialy) + ...
                                            T(i,j+1) * (k/h+kpartialy) + T(i,j-1) * (k/h-kpartialy) ...
                                        ) ...
                                    );
                            end 
                        end
                    
                    else 
                        if i == 1 || j == 1 || ... %Left and bottom sides
                             (i == n+1 && j <= n/2) || ... %Right and top sides
                             (i <= n/2 && j == n+1) || ... 
                             (i >= n/2 && j == n/2 ) || ... %Inward corner sides
                             (i == n/2 && j >= n/2 )
                             T(i,j) = 0;
                        elseif i > n/2 && j > n/2 
                             T(i,j) = missing;
                        else 
                             if nonLinearVersion == 0
                                 T(i,j) = h/(4*k)*...
                                    ( ...
                                        -2*h*f + ...
                                        ( ...
                                            T(i+1,j) * (k/h+kpartialx) + T(i-1,j) * (k/h-kpartialy) + ...
                                            T(i,j+1) * (k/h+kpartialy) + T(i,j-1) * (k/h-kpartialy) ...
                                        ) ...
                                    );
                             else 
                                T(i,j) = fsolve(@(T)( ...
                                         -( ...
                                         T(i+1,j) * (k/h+kpartialx) + T(i-1,j) * (k/h-kpartialy) + ...
                                         T(i,j+1) * (k/h+kpartialy) + T(i,j-1) * (k/h-kpartialy) ...
                                         -4*k/h*T ...
                                         ) ... 
                                         +2*h*T^4 - 2*h*f ...
                                    ), 5);
                             end
                        end
                    end
            end
            
        end
        
        noMissT = T(~ismissing(T)); noMissTemp = temp(~ismissing(T));
        residual = sqrt((sum(sum((noMissT-noMissTemp).^2)))/(n+1)^2);
        
    else 
            T = -T;
            break; 
    end
    
end 

%% Dumb way of getting # variables

k = 0; 
for i = 1:n+1 
    for j = 1:n+1
        if ~ismissing(T(i,j))
            k = k+1;
        end 
    end 
end 


%% Plot
%figure 
x = zeros(k,1);
y = zeros(k,1);
z = zeros(k,1);
idx = 1;
for i = 1:n+1 
    for j = 1:n+1
        if ~ismissing(T(i,j)) 
            x(idx) = i; y(idx) = j; z(idx) = T(i,j);
            idx = idx +1;
        end 
        
    end 
end 

tri = delaunay(x,y);
trisurf(tri,x,y,z)
        


