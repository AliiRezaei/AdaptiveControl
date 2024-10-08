clc
clear
close all

%% Simulation Time Vars

dt = 0.001;
SimTime = 5;
t = (0:dt:SimTime)';
nt = numel(t);

%% Single-Link Flexible Joint Robot System Params

Jm = 3.7e-3;
Jl = 9.3e-3;
m = 2.1e-1;
b = 3.1e-1 / 2;
k = 1.8e-1;
B = 4.6e-2;
k_tau = 8e-2;
g = 9.81;

%% Continuous Time Model

Ac = [   0,     1,     0,   0;
     -k/Jm, -B/Jm,  k/Jm,   0;
         0,     0,     0,   1;
      k/Jl,     0, -k/Jl,   0];

Bc = [0, k_tau/Jm, 0, 0]';

fc = @(x) [0, 0, 0, -m*g*b/Jl*sin(x)]';

% A_tilde = [-0.5, 0.5, 0, 0; -0.5, -0.5, 0, 0; 0, 0, -0.2, 0.3; 0, 0, -0.3, -0.2];

%% Discrete Time Model

A = eye(size(Ac)) + dt * Ac;
B = dt * Bc;
f = @(x) [0, 0, 0, -dt*m*g*b/Jl*sin(x)]';

%% Designing Params

M = zeros(4); M(4, 3) = sqrt(m*g*b/Jl);
N = zeros(4); N(3, 3) = sqrt(m*g*b/Jl) * dt;

H = [16.956, 0, 0, 0; 0, 17.82, -9.07, 0.138; 0, -9.07, 21.824, -1.526; 0, 0.138, -1.526, 0.4060];
Y = [0.8478, 0.0085, 0, 0; -0.4317, 0.8473, 0, 0; 0.2055, -0.0798, 0, 0; 0.0006, -0.0025, 0, 0];

Q = 0.05 * eye(4);
L = A - H \ Y;

A_tilde = [-0.5, 0.5, 0, 0; -0.5, -0.5, 0, 0; 0, 0, -0.2, 0.3; 0, 0, -0.3, -0.2];

%% Simulation

x = zeros(4, nt);
x(:, 1) = randn(4, 1);

x_bar = zeros(4, nt);
x_bar(:, 1) = randn(4, 1);

x_tilde = zeros(4, nt);
x_tilde(:, 1) = x(:, 1) - x_bar(:, 1);

theta = zeros(nt, 1);
theta_bar = zeros(nt, 1);
theta_tilde = zeros(nt, 1);

u = ones(nt, 1);
fault_time_min = 1;
fault_time_max = 3;
for kk = 1:nt - 1

    % u(kk) = 1;
    % u(kk) = B \ ((A_tilde - A) * x(:, kk) - f(x(3, kk))) / (1 - theta(kk));
    PI = - B * u(kk);

    if kk*dt > fault_time_min && kk*dt < fault_time_max
        theta(kk+1) = 2.5;
    end

    x(:, kk+1) = A * x(:, kk) + f(x(3, kk)) + B * u(kk) + PI * theta(kk);
    
    XI = 2 * (PI * PI' + Q)^(-1);
    
    x_tilde(:, kk+1) = L * x_tilde(:, kk) + PI * theta_tilde(kk) + (f(x(3, kk)) - f(x_bar(3, kk)));
    % x_tilde(:, kk+1) = x(:, kk) - x_bar(:, kk);
    x_bar(:, kk+1) = A * x_bar(:, kk) + f(x_bar(3, kk)) + B * u(kk) + PI * theta_bar(kk) + H \ Y * x_tilde(:, kk);
    
    theta_bar(kk+1) = theta_bar(kk) + PI' * XI * (x_tilde(:, kk+1) - L * x_tilde(:, kk) - (f(x(3, kk)) - f(x_bar(3, kk))));
    theta_tilde(kk+1) = theta(kk) - theta_bar(kk);

end

x = x';
x_bar = x_bar';
x_tilde = x_tilde';

figure
subplot(4, 1, 1)
plot(t, x(:, 1), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 1), '--', 'LineWidth', 2)
x_shade = [fault_time_min, fault_time_max, fault_time_max, fault_time_min];
y_shade = [min(x_bar(:, 1)), min(x_bar(:, 1)), max(x_bar(:, 1)), max(x_bar(:, 1))];
fill(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
legend('$True$ $State$', 'NLAO', 'Interpreter', 'latex')
ylabel('$x_{1k}$', 'Interpreter', 'latex')
title('$True$ $States$ $and$ $Observation$', 'Interpreter', 'latex')

subplot(4, 1, 2)
plot(t, x(:, 2), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 2), '--', 'LineWidth', 2)
y_shade = [min(x_bar(:, 2)), min(x_bar(:, 2)), max(x_bar(:, 2)), max(x_bar(:, 2))];
fill(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
ylabel('$x_{2k}$', 'Interpreter', 'latex')

subplot(4, 1, 3)
plot(t, x(:, 3), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 3), '--', 'LineWidth', 2)
y_shade = [min(x_bar(:, 3)), min(x_bar(:, 3)), max(x_bar(:, 3)), max(x_bar(:, 3))];
fill(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
ylabel('$x_{3k}$', 'Interpreter', 'latex')

subplot(4, 1, 4)
plot(t, x(:, 4), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 4), '--', 'LineWidth', 2)
y_shade = [min(x_bar(:, 4)), min(x_bar(:, 4)), max(x_bar(:, 4)), max(x_bar(:, 4))];
fill(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
xlabel('$k/\delta$ $(s)$', 'Interpreter', 'latex')
ylabel('$x_{4k}$', 'Interpreter', 'latex')

figure
plot(t, theta, 'LineWidth', 2)
hold on
plot(t, theta_bar, '--', 'LineWidth', 2)
xlabel('$k/\delta$ $(s)$', 'Interpreter', 'latex')
ylabel('$\theta$', 'Interpreter', 'latex')
title('$True$ $Fault$ $and$ $Estimation$', 'Interpreter', 'latex')
legend('$\theta$', '$\hat{\theta}$', 'Interpreter', 'latex')

figure
plot(t, x_tilde, 'LineWidth', 2)
xlabel('$k/\delta$ $(s)$', 'Interpreter', 'latex')
ylabel('$\tilde{x}_k$', 'Interpreter', 'latex')
title('$State$ $Estimation$ $Error$', 'Interpreter', 'latex')
legend('$\tilde{x}_1k$', '$\tilde{x}_2k$', '$\tilde{x}_3k$', '$\tilde{x}_4k$', 'Interpreter', 'latex')


