L = 27.68608264;
y = linspace (0,L,150); 
k = 0.056735954;
f = -(cos(k*y) - 1);
plot (f,y); 
xlim([0 3]);

title('Deflection of the Tree with Increasing Height');
xlabel('Deflection [m]'); 
ylabel('Height [m]'); 