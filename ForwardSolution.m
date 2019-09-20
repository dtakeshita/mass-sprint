%check the results via InverseProblem.m by solving the Newton's equation 
clear; close all;
load testdata.mat
ndat = length(savedata);
errYset = zeros(ndat,1);
errVset = zeros(ndat,1);
for nd = 1:10:ndat
    tmp = savedata{nd};
    tspan = [0 tmp.tfloor(end)];
    g = 9.8;
    X_CC0 = tmp.X_MTC0 - tmp.X_SEC0;
    F_SEC = tmp.M*(g-tmp.Afloor) - tmp.C*tmp.Vfloor;
    dX_SEC = F_SEC/tmp.K;
    X_SEC = tmp.X_SEC0 + dX_SEC;
    X_CC = tmp.Yfloor - X_SEC;
    %C = 1.98e3;
    ic = [tmp.X_MTC0; tmp.V0];
    opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
    [t,y] = ode45(@(t,y) focrceoscifcn(t,y,tmp.tfloor,X_CC,tmp.M,tmp.K,tmp.C,tmp.X_SEC0,g), tspan, ic, opts);
    %% forward solution
    omega_n = sqrt(tmp.K/tmp.M);
    omega_n2 = omega_n^2;
    Amp_CC = (max(X_CC) - min(X_CC))/2;
    omegaFwd = 2*pi/tmp.T;
    %Xfwd_CC = (Amp_CC)*(1 - cos(omegaFwd*tmp.tfloor)) + X_CC0;
    V0 = tmp.V0;
    Xfwd_CC = g*(1/omegaFwd^2 - 1/omega_n^2) ...
                - g*(1/omegaFwd^2 - 1/omega_n^2)*cos(omegaFwd*t)...
                + omegaFwd*V0*(1/omegaFwd^2 - 1/omega_n^2)*sin(omegaFwd*t)+ X_CC0;
    figure;
    plot(tmp.tfloor, X_CC)
    hold on
    plot(t, Xfwd_CC,'rx')    
    [tfwd, yfwd] = ode45(@(t,y) focrceoscifcn(t,y,tmp.tfloor,Xfwd_CC,tmp.M,tmp.K,tmp.C,tmp.X_SEC0,g), tspan, ic, opts);
    x_analy_nat1 = Amp_CC*omegaFwd^2/(omegaFwd^2 - omega_n2)*cos(omega_n*tmp.tfloor);
    x_analy_nat2 = tmp.V0/omega_n*sin(omega_n*tmp.tfloor);
    x_analy_sp = Amp_CC*omega_n2/(omega_n2 - omegaFwd^2)*cos(omegaFwd*tmp.tfloor);
    x_analy_const = tmp.X_MTC0 + tmp.M/tmp.K*g + Amp_CC + X_CC0;
    x_analy = x_analy_nat1 + x_analy_nat2 + x_analy_sp + x_analy_const;
    figure;
    plot(tfwd, yfwd(:,1),'b')
    hold on
    plot(t,y(:,1),'kx')
    plot(tmp.tfloor, x_analy,'rx')
    %plot(tmp.tfloor, -x_analy_sp + x_analy_const-X_CC0,'mo')

end




