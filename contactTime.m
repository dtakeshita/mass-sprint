clear; close all;
f_n = 2.5;%Hz
omega_n =  2*pi*f_n;
M = 1;
K = M*omega_n^2;
g = 9.8;
V0 = logspace(-2,1,100);
%% calculate contact time
ph = atan(-V0*omega_n/g)+pi;%
Tc = 2*ph/omega_n;
h0 = V0.^2/(2*g);
tTakeoffTh = 1/omega_n*(2*pi - acos((M*g-2*K*h0)./(M*g+2*K*h0)));
plot(V0,Tc)
hold on
plot(V0,tTakeoffTh,'x')
set(gca,'xscale','log')
