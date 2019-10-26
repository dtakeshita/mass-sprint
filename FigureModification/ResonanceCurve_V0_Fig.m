clear;
close all;
fig_dir = '/Users/dtakeshi/Dropbox/KyorinTextbook/FigureFIles/MassSpringWithCC';
%% Figure V0=0.1, T=0.35
fname = 'ResonanceCurve_V0.fig';
fh = open(fullfile(fig_dir,fname));
fwidth = 650;
fheight = 300;
fpos = [339   384   fwidth   fheight];
h_line = findobj(gca,'type','line');
xtick = [0.1 1];
xticklabel={'0.1','1'};
ytick = [1 10 100];
yticklabel={'1','10','100'};
xlim = [0.1 1];
ylim = [0.3 200];
set(h_line,'linewidth',2)
set(gca,'xscale','log','yscale','log','fontsize',18,...
    'xlim',xlim,'ylim',ylim,'xtick',xtick,'xticklabel',xticklabel,...
    'ytick',ytick','yticklabel',yticklabel)
c_max = 0.7;
n_clr = length(h_line);
c_set = linspace(0,c_max,n_clr)
for n = 1:n_clr
    c = c_set(n);
    set(h_line(n),'color',c*[1 1 1])
    
end
set(fh,'position',fpos)
