function [onBoundary] = isOnBoundary(i,j,n)
%Return boundary segment 
%   Return 0 if not not boundary. 1 if on boundary segments that are not
%   AB. 2 if on boundary segment AB
    if i == 1 || j == 1 || ... %Left and bottom sides
     (i == n+1 && j < n/2+1) || ... %Right and top sides
     (i < n/2+1 && j == n+1) || ...
     (j == n/2+1 && i > n/2 + 1 )
        onBoundary = 1;
    elseif (j > n/2 && i == n/2+1 ) %Inward corner sides
        onBoundary = 2;
    else
        onBoundary = 0;
    end
end

