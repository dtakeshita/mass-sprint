% calculate "contact time" in hopping of a simple mass-spring system 
%% ma = mg - F (downward is taken as positive)
clear;
close all;
m = 66.4;%kg
M = 520.38;%effective mass, see CalcEffectiveMass.m
%M = 200;
K = 1.657e5;%N/m - 165.7 N/mm
%C = 1.98e3;%N/(m/s) - 1.98 N/(mm/s)
C = 0;
omega_n = sqrt(K/M);
fn = omega_n/(2*pi);%Natural frequency in Hz

%omega_n = sqrt(K/M);
g = 9.80;

T_n = 2*pi/omega_n;
Tset = 0.45:0.01:1;
dt = 5e-4;

X_MTC0 = 0.5; %MTC length at rest
X_SEC0 = 0.3;
X_CC0 = X_MTC0 - X_SEC0;
V0set = [0];%initial downward velocity. must be zero or positive number 
tTakeoffSet = zeros(size(V0set));
h0 = V0set.^2/(2*g);
%tTakeoffTh = 1/omega*(2*pi - acos((M*g-2*k*h0)./(M*g+2*k*h0)));
ampRatio = zeros(size(Tset));
ndat = 1;
for n=1:length(V0set)
    ampRatio = zeros(size(Tset));
    tTakeoffSet = zeros(size(Tset));
    V0 = V0set(n);
    for nt = 1:length(Tset)
        T = Tset(nt);
        omega = 2*pi/T;
        t = 0:dt:1.1*T;
        a = 1;
        delta = -pi/4;
        Y1 = a*sin(omega*t + delta);

        Y = V0/omega*sin(omega*t) -g/omega^2 *cos(omega*t) + X_MTC0 + g/omega^2 ;
        V = V0*cos(omega*t) + g/omega*sin(omega*t);
        A = -V0*omega*sin(omega*t) + g*cos(omega*t);
        [tTakeoff, idxTakeoff] = findTcontact(t,V,T/2,-V0);
        Yshift = [Y(2:end) Y(end)];
        Yfloor = Y(1:idxTakeoff);
        Vfloor = V(1:idxTakeoff);
        Afloor = A(1:idxTakeoff);
        tfloor = t(1:idxTakeoff);
        F_SEC = M*(g-Afloor) - C*Vfloor;
        dX_SEC = F_SEC/K;
        X_SEC = X_SEC0 + dX_SEC;
        X_CC = Yfloor - X_SEC;
        dX_CC = X_CC - X_CC0;
        ampRatio(nt) = (max(Yfloor)-min(Yfloor))/(max(X_CC)-min(X_CC));
        tTakeoffSet(nt) = tTakeoff;
        %% save data and test it by solving ODE
        savedata{ndat} = struct('M',M,'K',K,'C',C,'X_MTC0',X_MTC0,'X_SEC0',X_SEC0,...
            'dt',dt,'T',T,'V0',V0,...
            'Yfloor',Yfloor,'Vfloor',Vfloor,'Afloor',Afloor,'tfloor',tfloor);
        ndat = ndat + 1;
        %% plot for the text book
        plot(tfloor, Yfloor-mean(Yfloor))
        hold on
        plot(tfloor, X_CC - mean(X_CC),'r')
        %% plot for checking purpose
%         subplot(6,1,1)
%         plot(tfloor,Yfloor)
%     %     hold on
%     %     plot(t,Y0*ones(size(t)),'r')
%     %     plot(tTakeoff, Y(idxTakeoff),'ko')
%         title('L_{MTC}')
%         subplot(6,1,2)
%         %plot(t,Y1)
%         plot(tfloor,Vfloor)
%     %     hold on
%     %     plot(t,-V0*ones(size(t)),'r')
%         title('V_{MTC}')
%         subplot(6,1,3)
%         plot(tfloor,Afloor)
%         title('Acceleration')
%         subplot(6,1,4)
%         plot(tfloor,F_SEC)
%         title('F_SEC')
%         subplot(6,1,5)
%         plot(tfloor,dX_SEC)
%         title('SEC displacement')
%         subplot(6,1,6)
%         plot(tfloor,dX_CC)
%         title('CC displacement')
    end
    save('testdata.mat', 'savedata');
    
    figure
    plot(tTakeoffSet, ampRatio)
    xlabel('Contact time (sec)')
    ylabel('A_{MTC}/A_{CC}')
    ttl_str = sprintf('Tn = %f, V0 = %f',T_n, V0);
    title(ttl_str)
end
% plot(h0,tTakeoffSet)
% hold on
% plot(h0, T*ones(size(h0)),'k')
% plot(h0, T/2*ones(size(h0)),'k')

