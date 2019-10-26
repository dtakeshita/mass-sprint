clear;
close all;
fig_dir = '/Users/dtakeshi/Dropbox/KyorinTextbook/FigureFIles/MassSpringWithCC';
%% Figure V0=0.1, T=0.35
fname = 'OptimalTc.fig';
fh = open(fullfile(fig_dir,fname));
fwidth = 650;
fheight = 300;
fpos = [339   384   fwidth   fheight];
h_line = findobj(gca,'type','line');
xtick = [0.01 0.1 1 10];
xticklabel={'0.01','0.1','1','10'};
ytick = 0.2:0.1:0.4;
%yticklabel={'0.1','0.2','1};
xlim = [0.01 10];
ylim = [0.15 0.4];
%set(h_line,'linewidth',2)
set(gca,'xscale','log','yscale','linear','fontsize',18,...
    'xlim',xlim,'ylim',ylim,'xtick',xtick,'xticklabel',xticklabel,...
    'ytick',ytick')
mkr_clr = 'k';
set(h_line(1),'MarkerSize',10,'MarkerEdgeColor',mkr_clr,...
    'MarkerFaceColor',mkr_clr)
set(h_line(2),'color','k','linestyle','-','linewidth',2)
set(fh,'position',fpos)
