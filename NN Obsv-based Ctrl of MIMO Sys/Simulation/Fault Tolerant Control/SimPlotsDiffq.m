clc
clear
close all
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaulttextinterpreter', 'latex')
figFontSize = 16; % figures font size
legFontSize = 13; % legends font size

%% System Spec

% validation end time
tf = 40; % seconds

% dimenssion
n = 2;
p = 1;


% control efort max/min
u_max =   1;
u_min = - 1;

%% Export Signals (temp)

% load data
out2 = load('SimResultsqIs2'); out2 = out2.out;
out15 = load('SimResultsqIs15'); out15 = out15.out;
out50 = load('SimResultsqIs50'); out50 = out50.out;

% time
t2   = out2.tout;
t15  = out15.tout;
t50  = out50.tout;


% control signal
uf      = out2.u_f.signals.values;
sat_u2  = out2.sat_u.signals.values;
sat_u15 = out15.sat_u.signals.values;
sat_u50 = out50.sat_u.signals.values;

%% Plots

% plot states & control signal
fig = figure;
theme(fig, 'light');
plot(t2, sat_u2); hold on;
plot(t15, sat_u15);
plot(t50, sat_u50);
plot(t2, uf, 'r');
yline(u_min, 'k--')
yline(u_max, 'k--')
ylim(1.15 * [u_min, u_max])
box on; grid on;
xlim([0 tf])
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('Control Signal', 'FontSize', figFontSize)
xlabel('$t$ [sec]', 'Interpreter', 'latex', 'FontSize', figFontSize)
legend('$u(t)$', '$u_f(t)$', '$u_{min/max}$', 'interpreter', 'latex', 'Location', 'northeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
set(gca, 'FontSize', figFontSize)

% % zoomed inset
% axes('Position',[0.24 0.218 0.21 0.09])
% plot(t, sat_u); hold on;
% plot(t, uf, 'r');
% yline(u_min, 'k--')
% yline(u_max, 'k--')
% box on; grid on;
% xlim([0 3])
% ylim(1.15 * [u_min, u_max])
% set(gca,'FontSize', 10)

% export high-quality
% set(gcf, 'Position', [100 100 700 850])
% exportgraphics(fig, 'Figures/XandU.pdf', 'ContentType', 'vector');
