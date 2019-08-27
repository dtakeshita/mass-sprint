%check the results via InverseProblem.m by solving the Newton's equation 
clear; close all;
load testdata.mat
ndat = length(savedata);
errYset = zeros(ndat,1);
errVset = zeros(ndat,1);
for nd = 1:ndat
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

    %ft = linspace(0,20,1000);
    %omega = 0.8*2*pi;
    %f = 0.1*sin(omega*ft);
    opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
    [t,y] = ode45(@(t,y) focrceoscifcn(t,y,tmp.tfloor,X_CC,tmp.M,tmp.K,tmp.C,tmp.X_SEC0,g), tspan, ic, opts);
    Ysolve = interp1(t,y(:,1),tmp.tfloor);
    Vsolve = interp1(t,y(:,2),tmp.tfloor);
    [errY, idx_errY] = max(abs(Ysolve-tmp.Yfloor)./tmp.Yfloor*100);
    [errV, idx_errV] = max(abs(Vsolve-tmp.Vfloor)./tmp.Vfloor*100);
%     [errY, idx_errY] = max(abs(Ysolve-tmp.Yfloor));
%     [errV, idx_errV] = max(abs(Vsolve-tmp.Vfloor));
    errYset(nd) = errY;
    errVset(nd) = errV;
    if errV > 5 %for small V (values close to 0) 
        figure;
        plot(t,y(:,1))
        hold on
        plot(tmp.tfloor,tmp.Yfloor,'rx')

        figure;
        plot(t,y(:,2))
        hold on
        plot(tmp.tfloor,tmp.Vfloor,'rx')
    end
end




