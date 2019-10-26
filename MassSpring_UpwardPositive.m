% calculate "contact time" in hopping of a simple mass-spring system 
% Upward is taken as the positive direction here
% v0 should be given as a positive number below
% ma = -mg - kx
% x = V0/omega_n sin(omega_n*t) - g/omega_n^2*cos(omega_n*t) + g/omega_n^2
% x(0) = 0, v0 = -|v0|, F0 = mg
clear;
close all;
M = 1;
g = 9.80;
dt = 1e-5;
Y0 = 0;
fn_set = [2.5 5];
fn_plot = [2.5];
V0_plot = [0.01 1.5199];%|V0| - given as a positive number
V0_plot = -V0_plot;
fh = figure;
fh_plot = figure;
for nf = 1:length(fn_set)
    fn = fn_set(nf);
    omega = 2*pi*fn;
    %omega = sqrt(K/M);
    K = M*omega^2;
    T = 2*pi/omega;
    t = 0:dt:1.1*T;
    V0set = logspace(-2,1,100);
    tTakeoffSet = zeros(size(V0set));
    h0 = V0set.^2/(2*g);
    tTakeoffTh = 1/omega*(2*pi - acos((M*g-2*K*h0)./(M*g+2*K*h0)));
    for n=1:length(V0set)
        V0 = -V0set(n);
        Y = V0/omega*sin(omega*t) + g/omega^2 *cos(omega*t) - g/omega^2;
        V = V0*cos(omega*t) - g/omega*sin(omega*t);
        Yshift = [Y(2:end) Y(end)];
        idxTakeoff = find(t > T/4 & Y < Y0 & Yshift > Y0, 1,'first');
%         Vshift = [V(2:end) V(end)];
%         idxTakeoff = find(t > T/4 & V < V0 & Vshift > V0, 1,'first');
%         if isempty(idxTakeoff)
%            idxTakeoff = find(t > T/4 & V == -V0, 1,'first'); 
%         end
        tTakeoff = t(idxTakeoff);
        tTakeoffSet(n) = tTakeoff;
        if any(abs(fn-fn_plot) < 1e-4) && any(abs(V0-V0_plot) < 1e-4)
            lw = 2;
            tplot = t(1:idxTakeoff);
            Yplot = Y(1:idxTakeoff);
            %tend = tplot(end);
            tend = 0.4
            tplot_end = ceil(tend/0.05)*0.05;
            xticks = 0:0.1:tplot_end;
%             if V0 > -0.011
%                 dy = 0.04;
%             else
%                 dy = 0.05;
%             end
            %ymin = floor(min(Yplot)/dy)*dy;
            dy = 0.05;
            ymin = -0.15;
            yticks = ymin:dy:0
            figure(fh_plot);
            h = plot(tplot,Yplot);
            hold on
            set(h,'linewidth',lw,'color','k','linestyle','--')
%             hold on
%             plot(tplot,Y0*ones(size(tplot)),'r')
%             plot(tTakeoff, Yplot(end),'ko')
            ttl_str = sprintf('fn=%g, V0=%g',fn, V0)
            title(ttl_str)
            set(gca,'fontsize',18,'xlim',[0 tplot_end ],'xtick',xticks,...
                'ylim',[ymin 0],'ytick',yticks)
        end
    end
    if nf==1
        cl = [0 0 0];
    else
        cl = [1 1 1]*0.6;
    end
    figure(fh)
    h = plot(V0set,tTakeoffSet);
    
    hold on
    set(h,'color',cl,'linewidth',lw)
    set(gca,'fontsize',18,'xscale','log',...
        'xticklabel',{'0.01';'0.1';'1';'10'})
    %plot(V0set, tTakeoffTh,'ro')
%     plot(V0set, T*ones(size(h0)),'k')
%     plot(V0set, T/2*ones(size(h0)),'k')
end


