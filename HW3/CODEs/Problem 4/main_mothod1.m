clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

problem.plant.Ap = [-6, 1, 0; -11, 0, 1; -6, 0, 0];
problem.plant.Bp = [0, 0, 1]';
problem.plant.Cp = [1, 0, 0];
problem.plant.kp = 1;
n  = size(problem.plant.Ap ,1);

problem.refModel.Am = [-12, 1, 0; -48, 0, 1; -64, 0, 0];
problem.refModel.Bm = [0, 0, 64]';
problem.refModel.Cm = [1, 0, 0];
problem.refModel.km = 64;
m  = size(problem.refModel.Am ,1);

problem.designParam.LAMBDA = - 1 * eye(n-1);
problem.designParam.l      =   1;
problem.designParam.gamma  =   10000;

%% Simulate System

dt = 0.01;                  % time step [seconds]
SimTime = 50;               % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)';    % time span
nt = numel(tSpan);
xp_0  = rand(n, 1);         % plant init cond
xm_0  = rand(m, 1);         % ref model init cond
w1_0 = zeros(n-1, 1);       % w1 init cond
w2_0 = zeros(n-1, 1);       % w2 init cond
phi_bar_0 = rand(2*n-1, 1); % estimated params init cond
e1_tilde_0 = zeros(m, 1);    % estimated params init cond
zeta_tilde_0 = zeros(m * (2*n-1), 1);    % estimated params init cond
InitCond = [xp_0; xm_0; w1_0; w2_0; phi_bar_0; e1_tilde_0; zeta_tilde_0]; % initial conditions
r = @(t) (1);               % step input
% r = @(t) sin(t);            % sine input
odeFunc = @(t, x) MRAC(t, x, r, problem); % ode function
[~, X] = ode45(odeFunc, tSpan, InitCond); % solve ode
nStates = size(X, 1);                     % number of states

% Unpack states :
xp        = X(:, 1:n);  
xm        = X(:, n+1:n+m);    
w1        = X(:, n+m+1:2*n+m-1);
w2        = X(:, 2*n+m:3*n+m-2);
phi_bar   = X(:, 3*n+m-1:5*n+m-3);

% Outputs :
yp = (problem.plant.Cp * xp')';
ym = (problem.refModel.Cm * xm')';

% Control signal :
w_bar = [w1, yp, w2];
u = sum((phi_bar .* w_bar), 2) + problem.plant.kp / problem.refModel.km * r(tSpan) .* ones(nt, 1);

%% Plots and Results

% plot yp and ym :
figure
plot(tSpan, yp, 'LineWidth', 2)
hold on
plot(tSpan, ym, '--', 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('y(t)', 'FontSize', 15)
title('$Plant$ $and$ $Reference$ $Model$ $Outputs$', 'FontSize', 15)
legend('$y_p(t)$','$y_m(t)$', 'interpreter', 'latex')

% plot estimated parameters :
figure
theta_1 = phi_bar(:, 1:n-1);
subplot(3, 1, 1)
plot(tSpan, theta_1, 'LineWidth', 2)
ylabel('$\hat{\theta}_1(t)$', 'FontSize', 15)
title('Estimated Parametes', 'FontSize', 15)

subplot(3, 1, 2)
theta_0 = phi_bar(:, n);
plot(tSpan, theta_0, 'LineWidth', 2)
ylabel('$\hat{\theta}_0(t)$', 'FontSize', 15)

subplot(3, 1, 3)
theta_2 = phi_bar(:, n+1:end);
plot(tSpan, theta_2, 'LineWidth', 2)
ylabel('$\hat{\theta}_2(t)$', 'FontSize', 15)
xlabel('t [sec]', 'FontSize', 15)

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
