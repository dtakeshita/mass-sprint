% simulation of d2Xdt2 = -k/mX - k/mAcc*cos(omega*t)
clear; close all;
g = 9.8;
%omega_drive = 1.3;
%f_drive = omega_drive/(2*pi);%driving frequency in Hz
f_drive = 120/60;%80/60, 120/60, 160/60, 200/60
C = 1.98e3;%N/(m/s) - 1.98 N/(mm/s)
%C = 0;
T = 1/f_drive;
dt = 1e-3;
M = 520.38;%effective mass, see CalcEffectiveMass.m
K = 1.657e5;%N/m - 165.7 N/mm

V0 = 0;
omegaFwd = 2*pi/T;
omega_n = sqrt(K/M);
tend = 50*T;
tspan = [0 tend];
tsim = 0:dt:tend;
%Amp_CC = 0.1;
%X0 = Amp_CC*K/M/(omega_n^2-omegaFwd^2);
X0 = 7;%In mm
%Amp_CC = M/K*X0*(omega_n^2-omegaFwd^2);
Amp_CC = sign(omega_n^2-omegaFwd^2)*X0*sqrt((1-(omegaFwd/omega_n)^2)^2+(omegaFwd*C/K)^2);
ic = [X0; V0];

Xfwd_CC = Amp_CC*cos(omegaFwd*tsim);
Xsp =1/(omega_n^2-omegaFwd^2)*K/M*Xfwd_CC;
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) forceoscifcn_simple(t,y,tsim,Xfwd_CC,M,K,C,g), tspan, ic, opts);
X_MTC = y(:,1);
X_CC = Xfwd_CC;
% plot the whold simulated time
plot(t,X_MTC,'b')
hold on
plot(tsim, X_CC)
%plot(tsim,Xsp)
legend({'MTC','CC','Xsp'})
%% cut only 1 cycle for plotting
Ncut = 20;
idx_use = t > Ncut*T;
X_use = X_MTC(idx_use);
t_use = t(idx_use);
d = diff(X_use);
idx_minima = find(d(1:end-1)  < 0 & d(2:end) > 0);
idx_plot = idx_minima(end-1)+1:idx_minima(end)+1;
t_plot = t_use(idx_plot);
X_plot = X_use(idx_plot);
idx_cc_plot = tsim >= t_plot(1) & tsim <=t_plot(end);
t_cc_plot = tsim(idx_cc_plot);
X_CC_plot = X_CC(idx_cc_plot);

% plot only 1 cycle
lw = 4;
t_x_norm = (t_plot-t_plot(1))/(t_plot(end)-t_plot(1))*100;
t_cc_norm = (t_cc_plot-t_cc_plot(1))/(t_cc_plot(end)-t_cc_plot(1))*100;
figure;
h_mtc = plot(t_x_norm, X_plot);
hold on
h_cc = plot(t_cc_norm, X_CC_plot);
set(h_mtc,'linewidth',lw,'color','k','linestyle','--')
set(h_cc,'linewidth',lw,'color','k','linestyle','-')
set(gca,'ylim',[-10 10],'xlim',[0 101],'xaxislocation','origin',...
    'ytick',-10:5:10,'xtick',0:50:100,'yticklabel',{},'xticklabel',{},...
    'linewidth',2)

if C == 0
    str_vis = 'No viscocity';
else
    str_vis = 'with viscocity';
end
str_ttl = sprintf('%s %.3gHz',str_vis,f_drive);
title(str_ttl)
