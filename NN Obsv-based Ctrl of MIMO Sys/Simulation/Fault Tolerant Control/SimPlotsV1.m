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
load SimResults4

% time
t  = out.tout; 
nt = numel(t);

% position
y     = out.y.signals.values;
x_hat = out.q_hat.signals.values;
x_des = out.q_des.signals.values;

% velocity
dx     = out.dq.signals.values;
dx_hat = out.dq_hat.signals.values;
dx_des = out.dq_des.signals.values;

% control signal
u      = out.u.signals.values;
uf     = out.u_f.signals.values;
sat_u  = out.sat_u.signals.values;

% fault spec
rho   = out.rho.signals.values;
delta = out.delta.signals.values;

% observer gain
eps = out.eps.signals.values;

% error signal
% e = out.e.signals.values;
e = reshape(out.e.signals.values, n, nt)';

% sliding surface
r_hat = out.r_hat.signals.values;

%% Plots

% plot states & control signal
fig = figure;
theme(fig, 'light');

subplot(3, 1, 1)
plot(t, y); hold on;
plot(t, x_hat, '--r');
plot(t, x_des, ':k');
box on; grid on;
xlim([0 tf])
ylabel('$y(t)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('System Output', 'FontSize', figFontSize)
legend('$y(t)$', '$\hat{x}_1(t)$', '$y_{d}(t)$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
set(gca, 'FontSize', figFontSize)

subplot(3, 1, 2)
plot(t, dx); hold on;
plot(t, dx_hat, '--r');
plot(t, dx_des, ':k');
box on; grid on;
xlim([0 tf])
ylim(1*[-2 1])
ylabel('$\dot{y}(t)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('System Output Derivative', 'FontSize', figFontSize)
legend('$\dot{y}(t)$', '$\hat{x}_2(t)$', '$\dot{y}_{d}(t)$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
set(gca, 'FontSize', figFontSize)

% zoomed inset
axes('Position',[0.24 0.44 0.11 0.09])
plot(t, dx); hold on;
plot(t, dx_hat, '--r');
plot(t, dx_des, ':k');
box on; grid on;
xlim([0 1])
set(gca,'FontSize', 10)

subplot(3, 1, 3)
plot(t, sat_u); hold on;
plot(t, uf, 'r');
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

% zoomed inset
axes('Position',[0.24 0.218 0.21 0.09])
plot(t, sat_u); hold on;
plot(t, uf, 'r');
yline(u_min, 'k--')
yline(u_max, 'k--')
box on; grid on;
xlim([0 3])
ylim(1.15 * [u_min, u_max])
set(gca,'FontSize', 10)

% export high-quality
set(gcf, 'Position', [100 100 700 850])
% exportgraphics(fig, 'Figures/XandU.pdf', 'ContentType', 'vector');

% plot error and sliding surface
fig = figure;
theme(fig, 'light');

subplot(3, 1, 1)
plot(t, r_hat);
box on; grid on;
xlim([0 tf])
ylim(1*[-2 1])
ylabel('$\hat{r}(t)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('Sliding Surface', 'FontSize', figFontSize)
set(gca, 'FontSize', figFontSize)

% zoomed inset
axes('Position',[0.22 0.74 0.11 0.09])
box on; grid on; hold on; 
plot(t, r_hat);
xlim([0 1])
set(gca,'FontSize', 10)

subplot(3, 1, 2)
plot(t, e(:, 1));
hold on; box on; grid on;
xlim([0 tf])
ylim(0.25*[-2 1])
ylabel('$y(t)-y_{d}(t)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('Output Tracking Error ($e$)', 'FontSize', figFontSize)
set(gca, 'FontSize', figFontSize)

% zoomed inset
axes('Position',[0.22 0.44 0.11 0.09])
plot(t, e(:, 1));
hold on; box on; grid on;
xlim([0 1])
set(gca,'FontSize', 10)

subplot(3, 1, 3)
plot(t, y - x_hat);
hold on; box on; grid on;
xlim([0 tf])
ylim(0.01*[-1 2])
ylabel('$y(t)-\hat{y}(t)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
xlabel('$t$ [sec]', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('Observation Error ($\tilde{y}$)', 'FontSize', figFontSize)
set(gca, 'FontSize', figFontSize)

% zoomed inset
axes('Position',[0.22 0.21 0.11 0.09])
plot(t, y - x_hat);
hold on; box on; grid on;
xlim([0 1])
% ylim(0.5*[-1 2])
set(gca,'FontSize', 10)

% export high-quality
set(gcf, 'Position', [100 100 700 850])
% exportgraphics(fig, 'Figures/Errors.pdf', 'ContentType', 'vector');
