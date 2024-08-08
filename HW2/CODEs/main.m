clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

problem.plant.Ap = [-2, 1; -2, 0];
problem.plant.Bp = [1, 2]';
problem.plant.Cp = [1, 0];
n  = size(problem.plant.Ap ,1);

problem.designParam.lambda =   10;
problem.designParam.LAMBDA = - 0.5;
problem.designParam.l      =   1;
problem.designParam.gamma  =   100;

%% Simulate System

dt = 0.01;                     % time step [seconds]
SimTime = 100;                 % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)';       % time span
x_0  = rand(n, 1);             % actual sys init cond
x1_hat_0 = rand(1, 1);         % estimated sys init cond
w1_hat_0 = rand(n-1, 1);       % w1 init cond
w2_hat_0 = rand(n-1, 1);       % w2 init cond
theta_hat_0 = rand(2*n, 1);    % estimated params init cond
InitCond = [x_0; x1_hat_0; w1_hat_0; w2_hat_0; theta_hat_0]; % initial conditions
% u = @(t) randn;                         % system input
% u = @(t) (1);                           % system input
% u = @(t) (sin(t) + sin(2*t) + sin(3*t));  % system input
u = @(t) (sin(t - pi/2) + sin(2*t + pi/3) + sin(3*t + pi/20) + sin(5*t + pi/2) + sin(7*t + pi/5));           % system input
odeFunc = @(t, x) NonMinimalObserver(t, x, u, problem); % ode function
[~, X] = ode45(odeFunc, tSpan, InitCond); % solve ode
nStates = size(X, 1);                     % number of states

% Unpack states :
x         = X(:, 1:n);    
x1_hat    = X(:, n+1);
w1_hat    = X(:, n+2:2*n);
w2_hat    = X(:, 2*n+1:3*n-1);
theta_hat = X(:, 3*n:end);

y = (problem.plant.Cp * x')';
y_hat = x1_hat;

%% Plots and Results

figure
plot(tSpan, y, 'LineWidth', 2)
hold on
plot(tSpan, y_hat, '--', 'LineWidth', 2)
title('$Measured$ $and$ $Estimated$ $Outputs$', 'FontSize', 15)
legend('$Measured$ $y$','$Estimated$ $y$', 'interpreter', 'latex')

% plot estimated parametes :
figure
c0 = theta_hat(:, 1);
subplot(4, 1, 1)
plot(tSpan, c0, 'LineWidth', 2)
ylabel('$\hat{c_0}(t)$', 'FontSize', 15)
title('Estimated Parametes', 'FontSize', 15)

subplot(4, 1, 2)
c = theta_hat(:, 2:n);
plot(tSpan, c, 'LineWidth', 2)
ylabel('$\hat{c}(t)$', 'FontSize', 15)

subplot(4, 1, 3)
d0 = theta_hat(:, n+1);
plot(tSpan, d0, 'LineWidth', 2)
ylabel('$\hat{d_0}(t)$', 'FontSize', 15)

subplot(4, 1, 4)
d = theta_hat(:, n+2:end);
plot(tSpan, d, 'LineWidth', 2)
xlabel('$t$ $[sec]$', 'FontSize', 15)
ylabel('$\hat{d}(t)$', 'FontSize', 15)

% Plot output error :
figure
plot(tSpan, y_hat - y, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('e(t)', 'FontSize', 15)
title('Output Error', 'FontSize', 15)

% Display estimated trandfer function :
s = tf('s');
PR = (c(end, :)) * inv(s * eye(n-1) - problem.designParam.LAMBDA) * [zeros(n-2, 1); problem.designParam.l] + c0(end); %#ok
QR = (d(end, :)) * inv(s * eye(n-1) - problem.designParam.LAMBDA) * [zeros(n-2, 1); problem.designParam.l] + d0(end); %#ok
W_s = minreal(PR / ((s + problem.designParam.lambda) - QR));
disp('Estimated Transfer Function is :')
W_s

