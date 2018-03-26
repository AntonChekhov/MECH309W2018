%% Initialize 

n = 8;

%Define matrices for parameters that need to be evaluated everywhere

T_ij = zeros(n+1, n+1);
%tau_ij evaluates T_ij term in the main PDE. 
%tau_im means tau_i-minus, means T_(i-1, j) term in main PDE. Similarly for
%j, and for p: tau_ip is T_(i+1, j) term in main PDE. 

%Initialize to same size
tau_ij = T_ij; tau_im = T_ij; tau_ip = T_ij; tau_jm = T_ij; tau_jp = T_ij;

%at this point I pay for starting my indices at 0 for all my derivations
for i = 1:n+1
    for j = 1:n+
        
        
        
        
    end
end