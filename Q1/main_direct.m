%% README
%This is attempt to directly construct equation system using the flatten
%and unflatten functions that relate grid positions to the "flattened"
%index in the temperature vector being solved for. 

%% CONFIG

n = 8;
h = 2/n;
maxIdx = flatten(n/2+1, n+1, n); %first two args give topmost righmost coord. Gives amount of vars
gradientBC = 0; %zero or one, depending on if one side should be insulated or not
nonLinearVersion = 0;  %zero or one, depending on if we solve nonlinear sys w/ fsolve
%Define matrices for parameters that need to be evaluated everywhere

%% System setup
A = zeros(2, maxIdx);



