%% README 
% file to test stuff

%% Flatten/unflatten
n = 8;
A = zeros(n+1, n+1);
v = ones(3, flatten(n/2+1, n+1, n)); 
k = 1;
for j = 1:n+1
    for i = 1:n+1
        A(i,j) = flatten(i,j,n);
        if ~ismissing(A(i,j))
            [tempi, tempj] = unflatten(A(i,j),n);
            v(1, k) = tempi;
            v(2, k) = tempj;
            v(3,k) = A(i,j);
            k = k+1;
        end 
        
    end 
end