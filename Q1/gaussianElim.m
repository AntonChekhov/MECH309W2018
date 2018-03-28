function [x] = gaussianElim(A, b)
%SOLVE Jh = -f(x_n) for Newton's method (using usual notation, Ax=b)
%   Using Gaussian elimination. 
n = length(A);
for j = 1:n-1
    if abs(A(j,j)/max(A(:,j))) < 0.01  %partial pivoting if needed
        [max_val , max_row] = max(A(:,j));
        temp = A(j, :); temp_b = b(j);
        A(j, :) = A(max_row, :); b(j) = b(max_row);
        A(max_row, :) = temp; 
    end

    for i = j+1:n   %standard gaussian elimination, forward
        mult = A(i,j)/A(j,j);
        for k = j:n
            A(i,k) = A(i,k) - mult*A(j,k); 
        end
        b(i) = b(i) - mult*b(j);
    end
end 

%backward substitution
x=zeros(1 , n);

for i = n:-1:1
    sum = 0;
    for j = n:-1:i+1
        sum = sum + x(j)*A(i,j);
    end 
    x(i) = (b(i) - sum)/A(i,i);
end 

x = x';
end





