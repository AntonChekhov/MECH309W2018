y = linspace(0,20,100);
%j = linspace(0,0,100);
k = 1;
for i = 0:0.001:10
    
    j(k) = besselFn(y(k),-1/3);
    k = k + 1;
end
plot (y,j)
axis([0,5,0,1])