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

%% Boundary
n = 8;
A = zeros(n+1, n+1);
v = ones(3, flatten(n/2+1, n+1, n)); 
k = 1;
for j = 1:n+1
    for i = 1:n+1
        A(i,j) = isOnBoundary(i,j,n);
        
    end 
end

%% System setup test
testVec = 1:maxIdx;
for k = 1:length(testVec)
    [i,j] = unflatten(k,n);
    testVec(k) = isOnBoundary(i,j,n);
    testVec = testVec ~= 0;
end

A2 = A(testVec, :);
A3 = sum(A2, :);