%% Initialize 
%n must be even for later stuff to compute. n/2 index badly handled
%otherwise
n = 8;
h = 1/n;
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
        f = k_ij(i*h, j*h);
        
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
            Tau_ij(i,j) = -1000000;
            Tau_ip(i,j) = -1000000;
            Tau_im(i,j) = -1000000;
            Tau_jp(i,j) = -1000000;
            Tau_jm(i,j) = -1000000;
            B_ij(i,j) = -1000000;
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
        
    end
end

