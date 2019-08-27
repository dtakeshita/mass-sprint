% calculate "contact time" in hopping of a simple mass-spring system 
clear;
close all;
M = 1;
K = 2;
g = 9.80;
omega = sqrt(K/M);
T = 2*pi/omega;
dt = 1e-3;
t = 0:dt:1.1*T;
Y0 = M*g/K;

V0set = -200:0.01:-0.01;
tTakeoffSet = zeros(size(V0set));
h0 = V0set.^2/(2*g);
tTakeoffTh = 1/omega*(2*pi - acos((M*g-2*K*h0)./(M*g+2*K*h0)));
for n=1:length(V0set)
    V0 = V0set(n);
    
    Y = V0/omega*sin(omega*t) + Y0 *cos(omega*t);
    V = V0*cos(omega*t) - g/omega*sin(omega*t);
    Yshift = [Y(2:end) Y(end)];
    idxTakeoff = find(t > T/4 & Y < Y0 & Yshift > Y0, 1,'first');
    tTakeoff = t(idxTakeoff);
    tTakeoffSet(n) = tTakeoff;
    
    subplot(2,1,1)
    plot(t,Y)
    hold on
    plot(t,Y0*ones(size(t)),'r')
    plot(tTakeoff, Y(idxTakeoff),'ko')
    subplot(2,1,2)
    plot(t,V)
    hold on
    plot(t,-V0*ones(size(t)),'r')
end
plot(h0,tTakeoffSet)
hold on
plot(h0, T*ones(size(h0)),'k')
plot(h0, T/2*ones(size(h0)),'k')

