clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

problem.plant.Ap = [3, 1; -2, 0];
problem.plant.Bp = [1, 1]';
problem.plant.Cp = [1, 0];
n  = size(problem.plant.Ap ,1);

problem.refModel.Am = - 1;
problem.refModel.Bm =   1;
problem.refModel.Cm =   1;
m  = size(problem.refModel.Am ,1);

problem.designParam.LAMBDA = - 1;
problem.designParam.l      =   1;
problem.designParam.gamma  =   100;

%% Simulate System

dt = 0.01;                     % time step [seconds]
SimTime = 20;                  % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)';       % time span
xp_0  = rand(n, 1);            % plant init cond
xm_0  = rand(m, 1);            % ref model init cond
w1_0 = rand(n-1, 1);           % w1 init cond
w2_0 = rand(n-1, 1);           % w2 init cond
theta_hat_0 = rand(2*n, 1);    % estimated params init cond
InitCond = [xp_0; xm_0; w1_0; w2_0; theta_hat_0]; % initial conditions
% r = @(t) (1);                             % step input
% r = @(t) sin(3*t);                             % sine input
r = @(t) (sin(t) + sin(2*t) + sin(3*t));  % sine input
odeFunc = @(t, x) MRAC(t, x, r, problem); % ode function
[~, X] = ode45(odeFunc, tSpan, InitCond); % solve ode
nStates = size(X, 1);                     % number of states

% Unpack states :
xp        = X(:, 1:n);  
xm        = X(:, n+1:n+m);    
w1        = X(:, n+m+1:2*n+m-1);
w2        = X(:, 2*n+m:3*n+m-2);
theta_hat = X(:, 3*n+m-1:end);

% Outputs :
yp = (problem.plant.Cp * xp')';
ym = (problem.refModel.Cm * xm')';

% Control signal :
w = [r(tSpan) .* ones(size(tSpan)), w1, yp, w2];
u = sum(theta_hat .* w, 2);

%% Plots and Results

% plot yp and ym :
figure
plot(tSpan, yp, 'LineWidth', 2)
hold on
plot(tSpan, ym, '--', 'LineWidth', 2)
title('$Plant$ $and$ $Reference$ $Model$ $Outputs$', 'FontSize', 15)
legend('$y_p(t)$','$y_m(t)$', 'interpreter', 'latex')

% plot estimated parameters :
figure
c0 = theta_hat(:, 1);
subplot(4, 1, 1)
plot(tSpan, c0, 'LineWidth', 2)
ylabel('$\hat{k}(t)$', 'FontSize', 15)
title('Estimated Parametes', 'FontSize', 15)

subplot(4, 1, 2)
c = theta_hat(:, 2:n);
plot(tSpan, c, 'LineWidth', 2)
ylabel('$\hat{\theta_1}(t)$', 'FontSize', 15)

subplot(4, 1, 3)
d0 = theta_hat(:, n+1);
plot(tSpan, d0, 'LineWidth', 2)
ylabel('$\hat{\theta_0}(t)$', 'FontSize', 15)

subplot(4, 1, 4)
d = theta_hat(:, n+2:end);
plot(tSpan, d, 'LineWidth', 2)
xlabel('$t$ $[sec]$', 'FontSize', 15)
ylabel('$\hat{\theta_2}(t)$', 'FontSize', 15)

% Plot output error :
figure
plot(tSpan, yp - ym, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('e(t)', 'FontSize', 15)
title('Output Error', 'FontSize', 15)

% Plot control signal :
figure
plot(tSpan, u, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('u(t)', 'FontSize', 15)
title('Control Signal', 'FontSize', 15)
