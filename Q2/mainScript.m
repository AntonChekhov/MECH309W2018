%main script to plot the deformed shapes

%values provided in the question
E=15000000000;
R=0.15;
S=pi*R^2;
I=(pi*R^4)/4;
rho=1000;
g=9.81;
q=rho*S*g/(E*I);

%Concentrated Load
LcConc = (pi^2/(4*q))^(1/3);
k = sqrt(rho*LcConc*S/(E*I));
y1 = linspace(0,LcConc,1000);
w1 = -(cos(k*y1) - 1);


%Distributed Load
kappasqrd=4*q/9;
root=bisection(1,5,0.00001,-1/3);
LcDist=(root^2 ./ kappasqrd)^(1/3);
y2 = linspace(0.01,LcDist,100); %y2 starts from 0.01 because the bessel function blows up at origin
for i = 1:100
    w2(i) = -1*RiemannSum(y2(i),LcDist,10) * 0.001;  %note that the function is scaled by 0.001
end

%plots
y2 = -y2 + LcDist;
plot(w1,y1,w2,y2);
axis([0,0.2,0,50])
xlabel('Deflection')
ylabel('Height')
title('Buckling Behaviour of a Tree')
legend ('Concentrated Load', 'Distributed Load')