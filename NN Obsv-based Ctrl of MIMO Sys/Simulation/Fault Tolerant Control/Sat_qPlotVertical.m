clc
clear
close all
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaulttextinterpreter', 'latex')
figFontSize = 28; % figures font size
legFontSize = 16; % legends font size
figTicksFontSize = 20;

%% Independent Var x

xLimit = [-8, 8];
nx     = 500;
x_min  = -2;
x_max  =  4;
x      = linspace(xLimit(1), xLimit(2), nx)';

%% Sat & Approximation

sat_q   = @(x, x_min, x_max, q) ((x_max+x_min)/2 + (x - (x_max+x_min)/2) ./ (1 + abs((x - (x_max+x_min)/2)/((x_max-x_min)/2)).^q).^(1/q));
dsat_q  = @(x, x_min, x_max, q) (1 + abs((x - (x_max+x_min)/2)/((x_max-x_min)/2)).^q).^(-(q+1)/q);
ddsat_q = @(x, x_min, x_max, p) ...
    -(p+1) * (x - (x_min + x_max)/2) .* abs(x - (x_min + x_max)/2).^(p-2) ...
    * ((x_max - x_min)/2)^(p+1) ./ ...
    ( ((x_max - x_min)/2)^p + abs(x - (x_min + x_max)/2).^p ).^((2*p+1)/p);

q  = [2, 5, 11.5];
nq = numel(q);

sat_x_val  = min(max(x, x_min), x_max);
dsat_x_val = 1 - min(1, max(x, x_min) - min(x, x_max));

sat_q_val   = zeros(nx, nq);
dsat_q_val  = zeros(nx, nq);
ddsat_q_val = zeros(nx, nq);

for iq = 1:nq
    sat_q_val(:, iq)   =   sat_q(x, x_min, x_max, q(iq));
    dsat_q_val(:, iq)  =  dsat_q(x, x_min, x_max, q(iq));
    ddsat_q_val(:, iq) = ddsat_q(x, x_min, x_max, q(iq));
end

%% Plots

fig = figure;
theme(fig, 'light');
tLayout = tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');
nexttile(1);

hold on; box on; grid on;
plot(x, sat_q_val(:, 1));
plot(x, sat_q_val(:, 2), '--');
plot(x, sat_q_val(:, 3), ':k');
ppp = plot(x, sat_x_val);
xline(0, 'k');
yline(0, 'k');
xlim(xLimit)
ylim([-5 5])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$\mathrm{sat}_q(x)$', 'Interpreter', 'latex')
title('Saturation Approximation')
legend(['$q=$', num2str(q(1))], ['$q=$', num2str(q(2))], ['$q=$', num2str(q(3))], '$q\to\infty$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 2, 'Orientation', 'horizontal', 'FontSize', legFontSize)

gca_instance = gca; 
set(gca_instance, 'FontSize', figFontSize); 
gca_instance.XAxis.FontSize = figTicksFontSize; 
gca_instance.XLabel.FontSize = figFontSize; 
gca_instance.YAxis.FontSize = figTicksFontSize; 
gca_instance.YLabel.FontSize = figFontSize;

nexttile(2); axis off;
nexttile(3); axis off;

% export high-quality
set(gcf, 'Position', [100 100 700 850])
exportgraphics(fig, 'Figures/Fig_Sim_Sat_q.pdf', 'ContentType', 'vector');


fig = figure;
theme(fig, 'light');
tLayout = tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');
nexttile(1);

hold on; box on; grid on;
plot(x, dsat_q_val(:, 1));
plot(x, dsat_q_val(:, 2), '--');
plot(x, dsat_q_val(:, 3), ':k');
plot(x, dsat_q(x, x_min, x_max, 1e7));
xline(0, 'k');
yline(0, 'k');
xlim(xLimit)
ylim(1.3*[-1/1.3 1])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$\frac{d}{dx}\mathrm{sat}_q(x)$', 'Interpreter', 'latex')
title('First Order Derivative')

gca_instance = gca; 
set(gca_instance, 'FontSize', figFontSize); 
gca_instance.XAxis.FontSize = figTicksFontSize; 
gca_instance.XLabel.FontSize = figFontSize; 
gca_instance.YAxis.FontSize = figTicksFontSize; 
gca_instance.YLabel.FontSize = figFontSize;

nexttile(2); axis off;
nexttile(3); axis off;

% export high-quality
set(gcf, 'Position', [100 100 700 850])
exportgraphics(fig, 'Figures/Fig_Sim_dSat_q.pdf', 'ContentType', 'vector');


fig = figure;
theme(fig, 'light');
tLayout = tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');
nexttile(1);

hold on; box on; grid on;
plot(x, ddsat_q_val(:, 1));
plot(x, ddsat_q_val(:, 2), '--');
plot(x, ddsat_q_val(:, 3), ':k');
plot([x_min x_min], [0 5], 'Color', ppp.Color);
plot([x_max x_max], [0 -5], 'Color', ppp.Color);
yline(0.01, 'Linewidth', 1.5, 'Color', ppp.Color)
yline(-0.01, 'Linewidth', 1.5, 'Color', ppp.Color)
xline(0, 'k');
yline(0, 'k');
xlim(xLimit)
ylim(1.3*[-1 1])
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$\frac{d^2}{dx^2}\mathrm{sat}_q(x)$', 'Interpreter', 'latex')
title('Second Order Derivative')

gca_instance = gca; 
set(gca_instance, 'FontSize', figFontSize); 
gca_instance.XAxis.FontSize = figTicksFontSize; 
gca_instance.XLabel.FontSize = figFontSize; 
gca_instance.YAxis.FontSize = figTicksFontSize; 
gca_instance.YLabel.FontSize = figFontSize;

nexttile(2); axis off;
nexttile(3); axis off;

% export high-quality
set(gcf, 'Position', [100 100 700 850])
exportgraphics(fig, 'Figures/Fig_Sim_ddSat_q.pdf', 'ContentType', 'vector');
