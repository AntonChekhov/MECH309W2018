y_nod = -1.5;
h = 0.001;
y = ones(1, 4/h); 
y(1) = y_nod;
k = 2;
for t = h:h:4
    y(k) = y(k-1)+h*yprime(y(k-1), t);
    k = k+1;
end

t = 0:h:4;
plot(t,y);