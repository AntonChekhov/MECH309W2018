function [val] = integrateSplineLength(t, k, h, N, xpt)
%INTEGRATESPLINELENGTH Summary of this function goes here
%   Detailed explanation goes here
sum = 0;
k_a = k(1:N); k_b = k(N+1:2*N); k_c = k(2*N+1:3*N);
h_a = h(1:N); h_b = k(N+1:2*N); h_c = k(2*N+1:3*N);

for i = 1:N
    for x = linspace(xpt(i), xpt(i+1), 10)
        sum = sum + integrand(t, k_b(i),k_c(i),h_b(i), h_c(i), x);
    end
end

val = sum;
end

function[val] = integrand(t, k_b,k_c,h_b, h_c, x)
    val = (k_b*t+h_b+2*(k_c*t+h_c)*x)*(k_b+2*k_c*x)/ ...
        sqrt( ...
            1+(k_b*t+h_b+2*(k_c*t+h_c)*x)^2 ...
            );
    
end
