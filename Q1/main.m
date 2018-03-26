%% Initialize 

n = 8;
h = 1/n;
%Define matrices for parameters that need to be evaluated everywhere

T_ij = zeros(n+1, n+1);
%tau_ij evaluates T_ij term in the main PDE. 
%tau_im means tau_i-minus, means T_(i-1, j) term in main PDE. Similarly for
%j, and for p: tau_ip is T_(i+1, j) term in main PDE. 

%Initialize to same size
tau_ij = T_ij; tau_im = T_ij; tau_ip = T_ij; tau_jm = T_ij; tau_jp = T_ij; 
b_ij = T_ij;

%Set up for square domain
for i = 1:n+1
    for j = 1:n+1
        [k, kpartialx, kpartialy] = k_ij((i-1)*h, (j-1)*h);
        f = k_ij(i*h, j*h);
        
        tau_ij(i,j) = -4*k/h;
        tau_ip(i,j) = k/h+kpartialx;
        tau_im(i,j) = k/h-kpartialx;
        tau_jp(i,j) = k/h+kpartialy;
        tau_jm(i,j) = k/h-kpartialy;
        b_ij(i,j) = 2*h*f;
    end
end

%Enforce BC's and put domain to L-shape

for i = 1:n+1
    for j = 1:n+1
        %Monkey values for points outside considered domain
        if i > n/2 && j > n/2
            tau_ij(i,j) = -1000000;
            tau_ip(i,j) = -1000000;
            tau_im(i,j) = -1000000;
            tau_jp(i,j) = -1000000;
            tau_jm(i,j) = -1000000;
            b_ij(i,j) = -1000000;
        end 
        
        %BC's
        if i == 0 || j == 0 ||
                
            
        end
        
    end
end

