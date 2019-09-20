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
Tset = 0.2:0.01:1;
Tplot = [0.2 0.35 1.0];
V0set = [0.01 0.1 1 10];%initial downward velocity. must be zero or positive number 
%V0set = 0:0.1:10;
V0plot = [0.1];
V0Plot_GainCurve = [0 0.01 0.1 1 10];
dt = 5e-4;

X_MTC0 = 0.5; %MTC length at rest
X_SEC0 = 0.3;
X_CC0 = X_MTC0 - X_SEC0;
tTakeoffSet = zeros(size(V0set));
h0 = V0set.^2/(2*g);
%tTakeoffTh = 1/omega*(2*pi - acos((M*g-2*k*h0)./(M*g+2*k*h0)));
ampRatio = zeros(size(Tset));
Tc_peak = zeros(size(V0set));
ndat = 1;
fh_ratio = figure;
for n=1:length(V0set)
    ampRatio = zeros(size(Tset));
    tTakeoffSet = zeros(size(Tset));
    V0 = V0set(n);
    for nt = 1:length(Tset)
        T = Tset(nt);
        omega = 2*pi/T;
        f = 1/T;
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
        if any(abs(T-Tplot) < 1e-4) && any(abs(V0-V0plot) < 1e-4)
            amp_MTC = (max(Yfloor)-min(Yfloor));
            amp_CC = (max(X_CC)-min(X_CC));
            amp_plot = max(amp_MTC, amp_CC);
            lw = 2;
            figure;
            Yplot = Yfloor-Yfloor(1);
            CCplot = X_CC - X_CC(1);
            h_MTC = plot(tfloor, Yplot);
            hold on
            h_CC = plot(tfloor, CCplot);
            set(h_MTC,'linewidth',lw,'color','k','linestyle','--')
            set(h_CC,'linewidth',lw,'color','k','linestyle','-')
            %[AX,H1,H2] = plotyy(tfloor, Yplot, tfloor, CCplot)
            str_ttl = sprintf('V0 =%g, T=%g',V0,T);
            title(str_ttl)
            yrange_MTC = [0 amp_plot];
            if T < T_n
                %yrange_CC = [X_CC0-amp_plot X_CC0];
                yrange_CC = [-amp_plot 0];
            else
                yrange_CC = [0 amp_plot];
            end
%             set(AX(1),'ylim',yrange_MTC)
%             set(AX(2),'ylim',yrange_CC)
        end
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
    [~,idx_max] = max(ampRatio);
    Tc_peak(n) = tTakeoffSet(idx_max);
    if any(abs(V0-V0Plot_GainCurve)<1e-3)
        figure(fh_ratio)
        plot(tTakeoffSet, ampRatio)
        %lh = plot(Tset, ampRatio)
        hold on
        xlabel('Contact time (sec)')
        ylabel('A_{MTC}/A_{CC}')
        ttl_str = sprintf('Tn = %f, V0 = %f',T_n, V0);
        title(ttl_str)
    end
    
end
V0_Tc = logspace(-2,1,100);
h0 = V0_Tc.^2/(2*g);
tTakeoffTh = 1/omega_n*(2*pi - acos((M*g-2*K*h0)./(M*g+2*K*h0)));
figure;
plot(V0_Tc,tTakeoffTh)
hold on
plot(V0set, Tc_peak,'o')
% plot(h0, T*ones(size(h0)),'k')
% plot(h0, T/2*ones(size(h0)),'k')
set(gca,'xscale','log')
xlabel('V0');ylabel('Optimal Tc')

