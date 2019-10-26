% calculate "contact time" in hopping of a simple mass-spring system 
% Downward is taken as the positive direction here
% ma = mg - kx
% x = V0/omega_n sin(omega_n*t) - g/omega_n^2*cos(omega_n*t) + g/omega_n^2
% x(0) = 0, v0 = +|v0|, F0 = mg
clear;
close all;
M = 650;
g = 9.80;
dt = 1e-5;
Y0 = 0;
fn_set = [2.5 5];
fn_plot = [2.5];
V0_plot = [0.01];
for nf = 1:length(fn_set)
    fn = fn_set(nf);
    omega = 2*pi*fn;
    %omega = sqrt(K/M);
    K = M*omega^2;
    T = 2*pi/omega;
    t = 0:dt:1.1*T;
    %V0set = 0.01:0.01:20;
    V0set = logspace(-2,1,100);
    tTakeoffSet = zeros(size(V0set));
    h0 = V0set.^2/(2*g);
    tTakeoffTh = 1/omega*(2*pi - acos((M*g-2*K*h0)./(M*g+2*K*h0)));
    for n=1:length(V0set)
        V0 = V0set(n);
        Y = V0/omega*sin(omega*t) - g/omega^2 *cos(omega*t) + g/omega^2;
        V = V0*cos(omega*t) + g/omega*sin(omega*t);
    %     Yshift = [Y(2:end) Y(end)];
    %     idxTakeoff = find(t > T/4 & Y < Y0 & Yshift > Y0, 1,'first');
        Vshift = [V(2:end) V(end)];
        idxTakeoff = find(t > T/4 & V < -V0 & Vshift > -V0, 1,'first');
        if isempty(idxTakeoff)
           idxTakeoff = find(t > T/4 & V == -V0, 1,'first'); 
        end
        tTakeoff = t(idxTakeoff);
        tTakeoffSet(n) = tTakeoff;
        if any(abs(fn-fn_set) < 1e-4) && any(abs(V0-V0_plot) < 1e-4)
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
    end
    plot(V0set,tTakeoffSet)
    hold on
    %plot(V0set, tTakeoffTh,'ro')
    plot(V0set, T*ones(size(h0)),'k')
    plot(V0set, T/2*ones(size(h0)),'k')
end
set(gca,'xscale','log')

