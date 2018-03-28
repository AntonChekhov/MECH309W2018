%% Initialize 
%n must be even for later stuff to compute. n/2 index badly handled
%otherwise
n = 8;
h = 2/n;
%Define matrices for parameters that need to be evaluated everywhere

T_ij = zeros(n+1, n+1);
%tau_ij evaluates T_ij term in the main PDE. 
%tau_im means tau_i-minus, means T_(i-1, j) term in main PDE. Similarly for
%j, and for p: tau_ip is T_(i+1, j) term in main PDE. 

%Initialize to same size
Tau_ij = T_ij; Tau_im = T_ij; Tau_ip = T_ij; Tau_jm = T_ij; Tau_jp = T_ij; 
B_ij = T_ij; F_ij = T_ij;

%Set up for square domain
for i = 1:n+1
    for j = 1:n+1
        [k, kpartialx, kpartialy] = k_ij((i-1)*h, (j-1)*h);
        f = f_ij(i*h, j*h);
        
        Tau_ij(i,j) = -4*k/h; %capitals denote matrices
        Tau_ip(i,j) = k/h+kpartialx;
        Tau_im(i,j) = k/h-kpartialx;
        Tau_jp(i,j) = k/h+kpartialy;
        Tau_jm(i,j) = k/h-kpartialy;
        B_ij(i,j) = 2*h*f;

    end
end

%Enforce BC's and put domain to L-shape

for i = 1:n+1
    for j = 1:n+1
        %Monkey values for points outside considered domain
        if i > n/2 && j > n/2
            Tau_ij(i,j) = missing;
            Tau_ip(i,j) = missing;
            Tau_im(i,j) = missing;
            Tau_jp(i,j) = missing;
            Tau_jm(i,j) = missing;
            B_ij(i,j) = missing;
            T_ij(i,j) = missing;
        end 
        
        %BC's - temperature zero all sides
        if i == 1 || j == 1 || ... %Left and bottom sides
             (i == n+1 && j <= n/2) || ... %Right and top sides
             (i <= n/2 && j == n+1) || ... 
             (i >= n/2 && j == n/2 ) || ... %Inward corner sides
             (i == n/2 && j >= n/2 )
                Tau_ij(i,j) = 1;
                Tau_ip(i,j) = 0;
                Tau_im(i,j) = 0;
                Tau_jp(i,j) = 0;
                Tau_jm(i,j) = 0;
                B_ij(i,j) = 0;
        end
        
        %Take away coefficients for boundaries where their corresponding
        %temperature term is undefined - e.g. Tau_im is undefined when i =
        %0
        if i == 1 
            Tau_im(i,j) = 0; %These go into the matrix. 
        end 
        
        if j == 1
            Tau_jm(i,j) = missing;
        end 
        
        if (i == n+1 && j <= n/2) || (i == n/2 && j >= n/2)
            Tau_ip(i,j) = 0;
        end 
        
        if (i <= n/2 && j == n+1) || (i > n/2 && j == n/2) 
            Tau_jp(i,j) = missing;
        end 
    
        
    end
end

%% Flatten created matrices to vectors 
%vectors are lowercase!!
%A(:) yields column concantenation of originals.
%Naming is unfortunate. LOWERCASE are flattened vectors. Capitalized are
%matrices. 

tau_ij = Tau_ij(:); 
tau_ip = Tau_ip(:);
tau_im = Tau_im(:);
tau_jp = Tau_jp(:);
tau_jm = Tau_jm(:);
b_ij = B_ij(:);
t_ij = B_ij(:);

clearvars -except  B_ij tau_ij tau_ip tau_im tau_jp tau_jm b_ij t_ij n T_ij Tau_ij Tau_jp Tau_jm Tau_im%clean up workspace cause it be messy

%Remove points outside boundaries from created vectors
tau_ij = tau_ij(~ismissing(tau_ij));
 
tau_ip = tau_ip(~ismissing(tau_ip)); tau_ip = tau_ip(1:(length(tau_ip)-1));
tau_im = tau_im(~ismissing(tau_im)); tau_im = tau_im(2:length(tau_im));

tau_jp = tau_jp(~ismissing(tau_jp));
tau_jm = tau_jm(~ismissing(tau_jm));

b_ij = b_ij(~ismissing(b_ij));
t_ij = t_ij(~ismissing(t_ij));
%% Construct and solve system
A = diag(tau_ij, 0) +diag(tau_im, -1) + diag(tau_ip,1) + ...
diag(tau_jm, -(n+1)); %+ diag(tau_jp, (n+1));
%Hotfix. The part corresponding to tau_jp that is corresponds to the right
%of that stupid BC is moved straight down to get it into the right
%equations. This is stupid and has to be changed for every change in matrix
%dimensions, should be better generalized later on. 
Tau_jp_hotfixed = diag(tau_jp, (n+1));
for i = 33:46
    Tau_jp_hotfixed(i+n/2+1, i+n+1) = Tau_jp_hotfixed(i, i+n+1);
    Tau_jp_hotfixed(i, i+n+1) = 0;
end 

A = A +Tau_jp_hotfixed;

A = -A; %To account for minus in the PDE LHS
b = b_ij;
%t = gaussianElim(A,b);
t = fsolve(@(x)(A*x-b), ones(56,1));
%% Reconstruct 2D domain for plotting

%i > n/2 && j > n/2 MISSING
k = 1; 
for i = 1:n+1 
    for j = 1:n+1
        if ~(i > n/2 && j > n/2)
            T_ij(i,j) = t(k);
            k = k+1;
        end 
    end 
end 
%% Delaunay n shit
figure 
k = 1; 
x = zeros(length(t),1);
y = zeros(length(t),1);
z = t';
for i = 1:n+1 
    for j = 1:n+1
        if ~ismissing(T_ij(i,j)) 
            x(k) = i; y(k) = j;
            k = k+1;
        end 
    end 
end 

tri = delaunay(x,y);
trisurf(tri,x,y,z)



