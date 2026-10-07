clc
clear
close all
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaulttextinterpreter', 'latex')
figFontSize = 28; % figures font size
legFontSize = 16; % legends font size
figTicksFontSize = 20;
figFontSizeZommed = 13;

%% System Spec

% validation end time
tf = 40; % seconds

% dimenssion
n   = 3;
rho = 2;
p   = 1;

% control efort max/min
u_max =   1.5;
u_min = - 1.5;

%% Export Signals

% load data
load SimResultsBenchmark

% time
t  = out.tout; 
nt = numel(t);

% position
y      = out.y.signals.values;
x_hat  = out.x_hat.signals.values;
y_des  = out.y_des.signals.values;
dy_des = out.dy_des.signals.values;
zeta   = reshape(out.zeta.signals.values, n, nt)';

% control signal
u      = out.u.signals.values;
uf     = out.u_f.signals.values;
sat_u  = out.sat_u.signals.values;

% fault spec
gamma = out.gamma.signals.values;
delta = out.delta.signals.values;

% observer gain
eps = out.eps.signals.values;

% error signal
e = reshape(out.e.signals.values, rho, nt)';
y_tilde = out.y_tilde.signals.values;

% sliding surface
r_hat = out.r_hat.signals.values;

%% Plots

% plot states & control signal
fig = figure;
theme(fig, 'light');

subplot(3, 1, 1)
plot(t, y); hold on;
plot(t, x_hat(:, 1), '--r');
plot(t, y_des, ':k');
box on; grid on;
xlim([0 tf])
ylim(1*[-3 1])
ylabel('$y(t)$', 'Interpreter', 'latex')
title('System Output')
legend('$y(t)$', '$\hat{x}_1(t)$', '$y_{d}(t)$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% zoomed inset
axes('Position',[0.24 0.74 0.11 0.09])
plot(t, y); 
hold on; box on; grid on;
plot(t, x_hat(:, 1), '--r');
plot(t, y_des, ':k');
xlim([0 1])
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSizeZommed);
set(gca_instance, 'XTick', [0 1]);
set(gca_instance, 'YTick', [str2double(gca_instance.YTickLabel{1}) str2double(gca_instance.YTickLabel{end})]);

subplot(3, 1, 2)
hold on;
plot(t, x_hat(:, 2), '--r');
plot(t, dy_des, ':k');
box on; grid on;
xlim([0 tf])
ylim(1*[-3 1])
ylabel('$\hat{x}_2(t)$', 'Interpreter', 'latex')
title('Observer Second State')
legend('$\hat{x}_2(t)$', '$\dot{y}_{d}(t)$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% zoomed inset
axes('Position',[0.24 0.44 0.11 0.09])
hold on; box on; grid on;
% plot(t, dx); hold on;
plot(t, x_hat(:, 2), '--r');
plot(t, dy_des, ':k');
xlim([0 1])
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSizeZommed);
set(gca_instance, 'XTick', [0 1]);
set(gca_instance, 'YTick', [str2double(gca_instance.YTickLabel{1}) str2double(gca_instance.YTickLabel{end})]);

subplot(3, 1, 3)
plot(t, zeta(:, 1)); hold on;
plot(t, zeta(:, 2), '--r');
plot(t, zeta(:, 3), ':k');
box on; grid on;
xlim([0 tf])
ylim(1*[-2 1])
ylabel('$\xi(t)$', 'Interpreter', 'latex')
xlabel('$t$ [sec]', 'Interpreter', 'latex')
title('System External/Internal Dynamics')
legend('$\xi_1(t)$', '$\xi_2(t)$', '$\eta(t)$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% export high-quality
set(gcf, 'Position', [100 100 700 850])
% exportgraphics(fig, 'Figures/Fig_Sim_States.pdf', 'ContentType', 'vector');

% plot error and sliding surface
fig = figure;
theme(fig, 'light');

subplot(3, 1, 1)
plot(t, r_hat);
box on; grid on;
xlim([0 tf])
ylim(0.5*[-2 1])
ylabel('$\hat{r}(t)$', 'Interpreter', 'latex')
title('Sliding Surface')
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% zoomed inset
axes('Position',[0.22 0.74 0.11 0.09])
hold on; box on; grid on;
plot(t, r_hat);
xlim([0 1])
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSizeZommed);
set(gca_instance, 'XTick', [0 1]);
set(gca_instance, 'YTick', [str2double(gca_instance.YTickLabel{1}) str2double(gca_instance.YTickLabel{end})]);

subplot(3, 1, 2)
plot(t, e(:, 1));
hold on; box on; grid on;
xlim([0 tf])
ylim(0.1*[-2 1])
ylabel('$y(t)-y_{d}(t)$', 'Interpreter', 'latex')
title('Output Tracking Error ($e$)')
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% zoomed inset
axes('Position',[0.22 0.44 0.11 0.09])
plot(t, e(:, 1));
hold on; box on; grid on;
xlim([0 1])
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSizeZommed);
set(gca_instance, 'XTick', [0 1]);
set(gca_instance, 'YTick', [str2double(gca_instance.YTickLabel{1}) str2double(gca_instance.YTickLabel{end})]);

subplot(3, 1, 3)
plot(t, y_tilde);
hold on; box on; grid on;
xlim([0 tf])
ylim(0.001*[-1 2])
ylabel('$y(t)-\hat{y}(t)$', 'Interpreter', 'latex')
xlabel('$t$ [sec]', 'Interpreter', 'latex')
title('Observation Error ($\tilde{y}$)')
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% zoomed inset
axes('Position',[0.22 0.16 0.11 0.09])
plot(t, y_tilde);
hold on; box on; grid on;
xlim([0 1])
% ylim(0.5*[-1 2])
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSizeZommed);
set(gca_instance, 'XTick', [0 1]);
set(gca_instance, 'YTick', [str2double(gca_instance.YTickLabel{1}) str2double(gca_instance.YTickLabel{end})]);

% export high-quality
set(gcf, 'Position', [100 100 700 850])
% exportgraphics(fig, 'Figures/Fig_Sim_Errors.pdf', 'ContentType', 'vector');

% plot control signal
fig = figure;
theme(fig, 'light');

subplot(3, 1, 1)
plot(t, sat_u); hold on;
plot(t, uf, 'r');
yline(u_min, 'k--')
yline(u_max, 'k--')
ylim(1.15 * [u_min, u_max])
box on; grid on;
xlim([0 tf])
ylabel('$u(t)$', 'Interpreter', 'latex')
title('Control Signal')
legend('$u(t)$', '$u_f(t)$', '$u_{min/max}$', 'interpreter', 'latex', 'Location', 'northeast')
legend('NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', legFontSize)
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

subplot(3, 1, 2)
plot(t, gamma);
box on; grid on;
xlim([0 tf])
ylabel('$\gamma(t)$', 'Interpreter', 'latex')
title('Loss of Effectivness Factor')
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

subplot(3, 1, 3)
plot(t, delta);
box on; grid on;
xlim([0 tf])
ylabel('$\delta(t)$', 'Interpreter', 'latex')
title('Additive Fault')
xlabel('$t$ [sec]', 'Interpreter', 'latex')
gca_instance = gca;
set(gca_instance, 'FontSize', figFontSize);
gca_instance.XAxis.FontSize = figTicksFontSize;
gca_instance.XLabel.FontSize = figFontSize;
gca_instance.YAxis.FontSize = figTicksFontSize;
gca_instance.YLabel.FontSize = figFontSize;

% export high-quality
set(gcf, 'Position', [100 100 700 850])
% exportgraphics(fig, 'Figures/Fig_Sim_ControlSignal.pdf', 'ContentType', 'vector');

