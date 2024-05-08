clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% The plant state-space as follows :
%   xp' = ap * xp + kp * u

problem.plant.ap = -1;      % plant feedback gain
problem.plant.kp =  1;      % plant feedforward gain

problem.adapt.am      = -1; % adaptive sys hurwitz param
problem.adapt.gamma   = 1; % adaptation rate
problem.adapt.modelID =  1; % set model ID for parametrization

%% Simulate System

dt = 0.1;                % time step [seconds]
SimTime = 30;            % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)'; % time span
xp_0  = rand;            % actual sys init cond
xp_hat_0 = rand;         % estimated sys init cond
ap_hat_0 = rand;         % estimated feedback gain init cond
kp_hat_0 = rand;         % estimated feedforward gain init cond
InitCond = [xp_0, xp_hat_0, ap_hat_0, kp_hat_0]'; % initial conditions
u = @(t) (1); % system input
odeFunc = @(t, x) AdaptIdent(t, x, u, problem); % ode function
[~, x] = ode45(odeFunc, tSpan, InitCond);       % solve ode

% Unpack states :
xp     = x(:, 1); % system measured state
xp_hat = x(:, 2); % system estimated state
ap_hat = x(:, 3); % estimated feedback gain
kp_hat = x(:, 4); % estimated feedforward gain

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
subplot(2, 1, 1)
plot(tSpan, ap_hat, 'LineWidth', 2)
ylabel('$\hat{a}_p(t)$', 'FontSize', 15)
title('Estimated Parameters', 'FontSize', 15)
subplot(2, 1, 2)
plot(tSpan, kp_hat, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$\hat{k}_p(t)$', 'FontSize', 15)

% Plot identification error :
figure
plot(tSpan, xp_hat - xp, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('e(t)', 'FontSize', 15)
title('Identification Error', 'FontSize', 15)
