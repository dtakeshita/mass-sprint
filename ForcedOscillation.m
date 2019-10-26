clear; close all;
% load testdata.mat
% nd = 1;
% tmp = savedata{nd};
g = 9.8;
%T = 0.4;
%T = tmp.T;
f_drive = 4;%driving frequency in Hz
T = 1/f_drive;
dt = 1e-3;
M = 520.38;%effective mass, see CalcEffectiveMass.m
%M = 200;
K = 1.657e5;%N/m - 165.7 N/mm
%C = 1.98e3;%N/(m/s) - 1.98 N/(mm/s)
C = 0;
%X_MTC0 = 0.5; %MTC length at rest
X_SEC0 = 0.3;
V0 = 0;
Amp_CC = 0.02;
omegaFwd = 2*pi/T;
omega_n = sqrt(K/M);
tend = 10*T;
tspan = [0 tend];
tsim = 0:dt:tend;
X_MTC0 = Amp_CC*omegaFwd^2/(omega_n^2-omegaFwd^2 )-M*g/K;
X_CC0 = X_MTC0 - X_SEC0;
ic = [X_MTC0; V0];
% X_CC0 = X_MTC0 - X_SEC0;
% F_SEC = M*(g-tmp.Afloor) - C*tmp.Vfloor;
% dX_SEC = F_SEC/K;
% X_SEC = X_SEC0 + dX_SEC;
% X_CC = tmp.Yfloor - X_SEC;
% Amp_CC = (max(X_CC) - min(X_CC))/2;
% Amp_CC_opt = (omega_n^2 - omegaFwd^2)/omegaFwd^2*M*g/K
Xfwd_CC = (-K/M)*Amp_CC*cos(omegaFwd*tsim) + X_CC0;
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) focrceoscifcn(t,y,tsim,Xfwd_CC,M,K,C,X_SEC0,g), tspan, ic, opts);
X_MTC = y(:,1)-mean(y(:,1));
X_CC = Xfwd_CC-X_CC0;
plot(t,X_MTC)
hold on
plot(tsim, X_CC)
legend({'MTC','CC'})