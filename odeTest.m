%check the results via InverseProblem.m by solving the Newton's equation 
clear; close all;
%load testdata.mat

tspan = [0 T];
ic = [0; 0];
ft = linspace(0,T,1000);
omega = 0.8*2*pi;
f = 0.1*sin(omega*ft);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) myODE(t,y,ft,f), tspan, ic, opts);
plot(t,y(:,1))
