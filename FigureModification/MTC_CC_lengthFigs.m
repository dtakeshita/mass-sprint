clear;
close all;
fig_dir = '/Users/dtakeshi/Dropbox/KyorinTextbook/FigureFIles/MassSpringWithCC';
%% Figure V0=0.1, T=0.35
fname = 'Length_MTC_CC_V0_1_T1_0.fig';
fh = open(fullfile(fig_dir,fname));
h_line = findobj(gca,'type','line');
t = get(h_line(1),'xdata');
x_lim = [0 t(end)];
y_lim = [-5 55];
ytick = 0:15:60;
xtick = 0:0.25:1.0;
ft_sz = 18;
set(gca,'fontsize',ft_sz,'xlim',x_lim,'ylim',y_lim,'ytick',ytick,'xtick',xtick)
%% Figure V0=0.1, T=0.2
% fname = 'Length_MTC_CC_V0_1_T0_2.fig';
% fh = open(fullfile(fig_dir,fname));
% h_line = findobj(gca,'type','line');
% t = get(h_line(1),'xdata');
% x_lim = [0 t(end)];
% ytick = -4:2:2;
% ft_sz = 18;
% set(gca,'fontsize',ft_sz,'xlim',x_lim,'ytick',ytick)

%% Figure V0=0.1, T=0.35
% fname = 'Length_MTC_CC_V0_1_T0_35.fig';
% fh = open(fullfile(fig_dir,fname));
% h_line = findobj(gca,'type','line');
% t = get(h_line(1),'xdata');
% x_lim = [0 t(end)];
% ytick = 0:2:6;
% xtick = 0:0.1:0.3;
% ft_sz = 18;
% set(gca,'fontsize',ft_sz,'xlim',x_lim,'ytick',ytick,'xtick',xtick)