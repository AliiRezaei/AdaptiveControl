clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% System Parameters

actualSys.A = [0, 1; -2, -3];
actualSys.B = [0, 1]';
n  = size(actualSys.A ,1);

adaptSys.Am = [0, 1; -1, -2];
adaptSys.Q  = eye(n);
adaptSys.P  = lyap(adaptSys.Am, adaptSys.Q);
adaptSys.gamma = 100;

%% Problem Definition

% The plant state-space as follows :
%   xp' = Ap * xp + Bp * u

problem.plant.ap    = -1;       % linear term gain
problem.plant.alpha = -2;       % nonlinear term gain
problem.plant.kp = 3;           % actual feedforward gain
problem.plant.f  = @(x) x^3;    % nonlinear function of xp [f(xp)]
problem.plant.g  = @(u) sin(u); % nonlinear function of  u [g(u)]

problem.adapt.am      = -2; % adaptive sys hurwitz param
problem.adapt.gamma   = 50; % adaptation rate
problem.adapt.modelID =  2; % set model ID for parametrization


%% Simulate System

dt = 0.1;                % time step [seconds]
SimTime = 60;            % maximum simulation time [seconds]
tSpan = (0:dt:SimTime)'; % time span
x_0  = rand(n, 1);       % actual sys init cond
x_hat_0 = rand(n, 1);    % estimated sys init cond
A_hat_0 = rand(n, n);    % estimated sys matrix init cond
B_hat_0 = rand(n, 1);    % estimated input matrix init cond
InitCond = [x_0; x_hat_0; A_hat_0(:); B_hat_0(:)]; % initial conditions
u = @(t) (sin(t) + sin(2*t) + sin(3*t));           % system input
odeFunc = @(t, x) AdaptIdentHighOrder(t, x, u, actualSys, adaptSys); % ode function
[~, X] = ode45(odeFunc, tSpan, InitCond); % solve ode
nStates = size(X, 1);                     % number of states

% Unpack states :
x     = X(:, 1:n);    
x_hat = X(:, n+1:2*n); 
A_hat = zeros(n, n, nStates);
c = 1;
for i = 1:n
    for j = 1:n
        A_hat(j, i, :) = X(:, 2*n+c);
        c = c + 1;
    end
end
B_hat = zeros(n, 1, nStates);
for k = 1:n
    B_hat(k, 1, :) = X(:, 2*n+c);
    c = c + 1;
end

%% Plots and Results

% plot x and x_hat :
figure
for i = 1:n
    plot(tSpan, x(:, i), 'LineWidth', 2)
    hold on
    plot(tSpan, x_hat(:, i), 'LineWidth', 2)
end
xlabel('$t$ $[sec]$', 'FontSize', 15)
ylabel('$x(t)$', 'FontSize', 15)
title('Measured and Estimated States', 'FontSize', 15)
legend('$x_1(t)$','$\hat{x}_1(t)$', ...
    '$x_2(t)$','$\hat{x}_2(t)$', ...
    'interpreter', 'latex')

% plot estimated parametes :
figure
subplot(2,1,1)
for i = 1:n
    for j = 1:n
        plot(tSpan, squeeze(A_hat(i, j, :)), 'LineWidth', 2)
        hold on
    end
end
title('$\hat{A}(t)$', 'FontSize', 15)
legend('$\hat{A}_{11}$', '$\hat{A}_{12}$', ...
        '$\hat{A}_{21}$', '$\hat{A}_{22}$', ...
        'interpreter', 'latex')

subplot(2,1,2)
for i = 1:n
    plot(tSpan, squeeze(B_hat(i, 1, :)), 'LineWidth', 2)
    hold on
end
xlabel('$t$ $[sec]$', 'FontSize', 15)
title('$\hat{B}(t)$', 'FontSize', 15)
legend('$\hat{B}_{11}$', '$\hat{B}_{21}$', ...
        'interpreter', 'latex')

% Display final estimation :
disp('Estimated A Matrix is :')
disp(A_hat(:, : , end))
disp('Estimated B Matrix is :')
disp(B_hat(:, 1 , end))