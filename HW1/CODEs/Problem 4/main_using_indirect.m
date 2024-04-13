clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% The plant state-space as follows :
%   xp' = ap * xp + kp * u
% The reference model state-space as follows :
%   xm' = am * xm + km * r

problem.plant.ap = 2; % plant feedback gain
problem.plant.kp = 2; % plant feedforward gain

problem.refModel.am    =  -3; % ref model feedback gain
problem.refModel.km    =   3; % ref model feedforward gain
problem.refModel.gamma =  50; % adaptation rate

%% Simulate System

dt = 0.1;                % time step [seconds]
SimTime = 30;            % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)'; % time span
xp_0  = rand;            % actual sys init cond
xp_hat_0  = rand;        % estimated sys init cond
xm_0 = rand;             % ref model init cond
ap_hat_0 = rand;         % estimated feedback gain init cond
kp_hat_0 = 3 + rand;     % estimated feedforward gain init cond
InitCond = [xp_0, xp_hat_0, xm_0, ap_hat_0, kp_hat_0]'; % initial conditions
r = @(t) sin(t);                                  % reference signal
odeFunc = @(t, x) IndirectMRAC(t, x, r, problem); % ode function
[~, x] = ode45(odeFunc, tSpan, InitCond);         % solve ode

% Unpack states :
xp     = x(:, 1); % plant measured state
xp_hat = x(:, 2); % plant estimated state
xm     = x(:, 3); % ref model state
ap_hat = x(:, 4); % plant estimated feedback gain
kp_hat = x(:, 5); % plant estimated feedforward gain

%% Plots and Results

% Plot xp and xm :
figure
plot(tSpan, xp, 'LineWidth', 2)
hold on
plot(tSpan, xp_hat, 'LineWidth', 2)
plot(tSpan, xm, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('x(t)', 'FontSize', 15)
title('Reference Model and Plant States (Estimated and Measured)', 'FontSize', 15)
legend('$x_p(t)$','$\hat{x}_p(t)$','$x_m(t)$','interpreter', 'latex')

% Plot estimated parametes :
figure
subplot(2, 1, 1)
plot(tSpan, ap_hat, 'LineWidth', 2)
ylabel('$\hat{a}_p(t)$', 'FontSize', 15)
title('Estimated Plant Parameters', 'FontSize', 15)
subplot(2, 1, 2)
plot(tSpan, kp_hat, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$\hat{k}_p(t)$', 'FontSize', 15)

% Plot identification error :
figure
plot(tSpan, xp_hat - xp, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$e_i(t)$', 'FontSize', 15)
title('Identification Error', 'FontSize', 15)

% Plot tracking error :
figure
plot(tSpan, xp - xm, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$e_t(t)$', 'FontSize', 15)
title('Tracking Error', 'FontSize', 15)

% Plot control signal :
figure
u = (problem.refModel.am - ap_hat) ./ (kp_hat) .* xp + ...
    (problem.refModel.km) ./ (kp_hat) .* r(tSpan);
plot(tSpan, u, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('u(t)', 'FontSize', 15)
title('Control Signal', 'FontSize', 15)
