clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% The plant state-space as follows :
%   xp' = ap * xp + alpha * f(xp) + kp * g(u)

problem.plant.ap    = -1;       % linear term gain
problem.plant.alpha = -2;       % nonlinear term gain
problem.plant.kp = 3;           % actual feedforward gain
problem.plant.f  = @(x) x^3;    % nonlinear function of xp [f(xp)]
problem.plant.g  = @(u) sin(u); % nonlinear function of  u [g(u)]

problem.adapt.am      = -1; % adaptive sys hurwitz param
problem.adapt.gamma   = 10; % adaptation rate
problem.adapt.modelID =  1; % set model ID for parametrization

%% Simulate System

dt = 0.1;                % time step [seconds]
SimTime = 30;            % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)'; % time span
xp_0  = rand;            % actual sys init cond
xp_hat_0 = rand;         % estimated sys init cond
ap_hat_0 = rand;         % estimated linear term gain init cond
alpha_hat_0 = rand;      % estimated nonlinear term gain init cond
kp_hat_0 = rand;         % estimated feedforward gain init cond
InitCond = [xp_0, xp_hat_0, ap_hat_0, alpha_hat_0, kp_hat_0]'; % initial conditions
u = @(t) sin(t) + 0.2 * sin(2 * t);          % system input
odeFunc = @(t, x) AdaptIdentNonLin(t, x, u); % ode function
[~, x] = ode45(odeFunc, tSpan, InitCond);    % solve ode

% Unpack states :
xp        = x(:, 1); % system measured state
xp_hat    = x(:, 2); % system estimated state
ap_hat    = x(:, 3); % estimated linear term gain
alpha_hat = x(:, 4); % estimated nonlinear term gain
kp_hat    = x(:, 5); % estimated feedforward gain

%% Plots and Results

% Plot x and xp_hat :
figure
plot(tSpan, xp, 'LineWidth', 2)
hold on
plot(tSpan, xp_hat, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('x(t)', 'FontSize', 15)
title('Measured and Estimated States', 'FontSize', 15)
legend('$x_p(t)$','$\hat{x}_p(t)$','interpreter', 'latex')

% Plot estimated parametes :
figure
subplot(3, 1, 1)
plot(tSpan, ap_hat, 'LineWidth', 2)
ylabel('$\hat{a}_p(t)$', 'FontSize', 15)
title('Estimated Parameters', 'FontSize', 15)

subplot(3, 1, 2)
plot(tSpan, alpha_hat, 'LineWidth', 2)
ylabel('$\hat{\alpha}(t)$', 'FontSize', 15)

subplot(3, 1, 3)
plot(tSpan, kp_hat, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$\hat{k}_p(t)$', 'FontSize', 15)
