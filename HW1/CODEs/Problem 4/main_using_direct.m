clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% The plant state-space as follows :
%   xp' = ap * xp + kp * u
% The reference model state-space as follows :
%   xm' = am * xm + km * r

problem.plant.ap = 2; % actual sys feedback gain
problem.plant.kp = 2; % actual sys feedforward gain

problem.refModel.am    =  -3; % ref model feedback gain
problem.refModel.km    =   3; % ref model feedforward gain
problem.refModel.gamma = 100; % adaptation rate

%% Simulate System

dt = 0.1;                % time step [seconds]
SimTime = 30;            % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)'; % time span
xp_0  = rand;            % actual sys init cond
xm_0 = rand;             % estimated sys init cond
theta_hat_0 = rand;      % estimated feedback gain init cond
k_hat_0 = rand;          % estimated feedforward gain init cond
InitCond = [xp_0, xm_0, theta_hat_0, k_hat_0]'; % initial conditions
r = @(t) sin(t) + sin(2*t);                                % reference signal
odeFunc = @(t, x) DirectMRAC(t, x, r, problem);          % ode function
[~, x] = ode45(odeFunc, tSpan, InitCond);       % solve ode

% Unpack states :
xp        = x(:, 1); % plant measured state
xm        = x(:, 2); % ref model state
theta_hat = x(:, 3); % controller param
k_hat     = x(:, 4); % controller param

%% Plots and Results

% Plot xp and xm :
figure
plot(tSpan, xp, 'LineWidth', 2)
hold on
plot(tSpan, xm, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('x(t)', 'FontSize', 15)
title('Reference Model and Plant States', 'FontSize', 15)
legend('$x_p(t)$','$x_m(t)$','interpreter', 'latex')

% Plot estimated parametes :
figure
subplot(2, 1, 1)
plot(tSpan, theta_hat, 'LineWidth', 2)
ylabel('$\hat{\theta}(t)$', 'FontSize', 15)
title('Estimated Controller Parameters', 'FontSize', 15)
subplot(2, 1, 2)
plot(tSpan, k_hat, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$\hat{k}(t)$', 'FontSize', 15)
