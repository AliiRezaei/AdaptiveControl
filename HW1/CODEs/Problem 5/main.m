clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% The plant state-space as follows :
%   xp' = Ap * xp + Bp * u
% The reference model state-space as follows :
%   xm' = Am * xm + Bm * r

problem.plant.Ap = [0, 1; -4, 4];
problem.plant.Bp = [0, 1]';
n  = size(problem.plant.Ap ,1);

problem.refModel.Am = [0, 1; -1, -1];
problem.refModel.Bm = [0, 1]';
Q  = eye(n);
problem.refModel.P  = lyap(problem.refModel.Am, Q);
problem.refModel.gamma = 1;

%% Simulate System

dt = 0.1;                        % time step [seconds]
SimTime = 60;                    % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)';         % time span
xp_0  = rand(n, 1);              % plant init cond
xm_0  = xp_0 + 0.1 * rand(n, 1); % ref model init cond
theta_hat_0 = rand(n, 1);
InitCond = [xp_0; xm_0; theta_hat_0];              % initial conditions
r = @(t) (sin(t) + sin(2*t) + sin(3*t));           % system input
odeFunc = @(t, x) HighOrderMRAC(t, x, r, problem); % ode function
[~, X] = ode45(odeFunc, tSpan, InitCond);          % solve ode
nStates = size(X, 1);                              % number of states

% Unpack states :
x_p       = X(:, 1:n);    
x_m       = X(:, n+1:2*n); 
theta_hat = X(:, 2*n+1:3*n);

%% Plots and Results

% plot x and x_hat :
figure
plot(tSpan, x_p, 'LineWidth', 2)
hold on
plot(tSpan, x_m, 'LineWidth', 2)
xlabel('$t$ $[sec]$', 'FontSize', 15)
ylabel('$x(t)$', 'FontSize', 15)
title('Reference Model and Plant States', 'FontSize', 15)
legend('$x_{p1}(t)$','$x_{p2}(t)$', ...
    '$x_{m1}(t)$','$x_{m1}(t)$', ...
    'interpreter', 'latex')

% plot estimated parametes :
figure
plot(tSpan, theta_hat, 'LineWidth', 2)
ylabel('$\hat{\theta}(t)$', 'FontSize', 15)
title('Estimated Controller Parameters', 'FontSize', 15)
legend('$\hat{\theta}_{1}$', '$\hat{\theta}_{2}$', 'interpreter', 'latex')

% Plot tracking error :
figure
plot(tSpan, x_p - x_m, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('$e_t(t)$', 'FontSize', 15)
title('Tracking Error', 'FontSize', 15)

% Plot control signal :
figure
u = zeros(nStates, n);
for k = 1:nStates
    u(k, :) = theta_hat(k, :) * x_p(k, :)' + r(tSpan(k));
end
plot(tSpan, u, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('u(t)', 'FontSize', 15)
title('Control Signal', 'FontSize', 15)
