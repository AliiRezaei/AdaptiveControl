clc
clear
close all
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaulttextinterpreter', 'latex')
figFontSize = 16; % figures font size
legFontSize = 13; % legends font size

%% Independent Var x

xLimit = [-8, 8];
nx     = 500;
x_min  = -2;
x_max  =  4;
x      = linspace(xLimit(1), xLimit(2), nx)';

%% Sat & Approximation

sat_q  = @(x, x_min, x_max, q) ((x_max+x_min)/2 + (x - (x_max+x_min)/2) ./ (1 + abs((x - (x_max+x_min)/2)/((x_max-x_min)/2)).^q).^(1/q));
dsat_q = @(x, x_min, x_max, q) (1 + abs((x - (x_max+x_min)/2)/((x_max-x_min)/2)).^q).^(-(q+1)/q);

q  = [2, 5, 11.5];
nq = numel(q);

sat_x_val  = min(max(x, x_min), x_max);
dsat_x_val = 1 - min(1, max(x, x_min) - min(x, x_max));

sat_q_val  = zeros(nx, nq);
dsat_q_val = zeros(nx, nq);

for iq = 1:nq
    sat_q_val(:, iq)  =  sat_q(x, x_min, x_max, q(iq));
    dsat_q_val(:, iq) = dsat_q(x, x_min, x_max, q(iq));
end

%% Plots

fig = figure;
theme(fig, 'light');
colors = hsv(nq+1);
subplot(2, 1, 1)
hold on; box on; grid on;
plot(x, sat_q_val(:, 1));
plot(x, sat_q_val(:, 2), '--');
plot(x, sat_q_val(:, 3), ':k');
plot(x, sat_x_val);
xline(0, 'k');
yline(0, 'k');
xlim(xLimit)
ylim([-5 5])
ylabel('$\mathrm{sat}_q(x)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('Saturation Approximation', 'FontSize', figFontSize)
legend(['$q=$', num2str(q(1))], ['$q=$', num2str(q(2))], ['$q=$', num2str(q(3))], '$q\to\infty$', 'interpreter', 'latex', 'Location', 'southeast')
legend('NumColumns', 2, 'Orientation', 'horizontal', 'FontSize', legFontSize)
set(gca, 'FontSize', figFontSize)

subplot(2, 1, 2)
hold on; box on; grid on;
plot(x, dsat_q_val(:, 1));
plot(x, dsat_q_val(:, 2), '--');
plot(x, dsat_q_val(:, 3), ':k');
plot(x, dsat_q(x, x_min, x_max, 1e7));
xline(0, 'k');
yline(0, 'k');
xlim(xLimit)
ylim(1.3*[-1/1.3 1])
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', figFontSize)
ylabel('$\frac{d}{dx}\mathrm{sat}_q(x)$', 'Interpreter', 'latex', 'FontSize', figFontSize)
title('Saturation Approximation Derivative', 'FontSize', figFontSize)
set(gca, 'FontSize', figFontSize)

% export high-quality
set(gcf, 'Position', [100 100 700 850])
exportgraphics(fig, 'Figures/Sat_q.pdf', 'ContentType', 'vector');
