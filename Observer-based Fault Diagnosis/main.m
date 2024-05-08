clc
clear
close all

%% Simulation Time Vars

dt = 0.001;
SimTime = 5;
t = (0:dt:SimTime)';
nt = numel(t);

%% Single-Link Flexible Joint Robot System Params

A = [1, 1e-3, 0, 0; -0.05, 1, 0.05, 0; 0, 0, 1, 1e-3; 0.02, 0, -0.02, 1];
B = [0, 0.0216, 0, 0]';
f = @(x) [0, 0, 0, -0.0333 * sin(x)]';

% A_tilde = diag(linspace(-0.8, 0.8, 4));
% A_tilde = 0.1 * eye(4);
A_tilde = [-0.5, 0.5, 0, 0; -0.5, -0.5, 0, 0; 0, 0, -0.2, 0.3; 0, 0, -0.3, -0.2];

%% Designing Params

M = zeros(4); M(4, 3) = 0.06;
N = zeros(4); N(3, 3) = 0.06;

H = [16.956, 0, 0, 0; 0, 17.82, -9.07, 0.138; 0, -9.07, 21.824, -1.526; 0, 0.138, -1.526, 0.4060];
Y = [0.8478, 0.0085, 0, 0; -0.4317, 0.8473, 0, 0; 0.2055, -0.0798, 0, 0; 0.0006, -0.0025, 0, 0];

Q = 0.5 * eye(4);
L = A - H \ Y;

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

for k = 1:nt - 1

    % u(k) = 1;
    % u(k) = B \ ((A_tilde - A) * x(:, k) - f(x(3, k))) / (1 - theta(k));
    PI = - B * u(k);

    if k*dt > 2
        theta(k+1) = 1.5;
    end

    x(:, k+1) = A * x(:, k) + f(x(3, k)) + B * u(k) + PI * theta(k);
    
    XI = 2 * (PI * PI' + Q)^(-1);
    
    x_tilde(:, k+1) = L * x_tilde(:, k) + PI * theta_tilde(k) + (f(x(3, k)) - f(x_bar(3, k)));
    
    x_bar(:, k+1) = A * x_bar(:, k) + f(x_bar(3, k)) + B * u(k) + PI * theta_bar(k) + H \ Y * x_tilde(:, k);
    
    theta_bar(k+1) = theta_bar(k) + PI' * XI * (x_tilde(:, k+1) - L * x_tilde(:, k) - (f(x(3, k)) - f(x_bar(3, k))));
    theta_tilde(k+1) = theta(k) - theta_bar(k);

end

x = x';
x_bar = x_bar';
x_tilde = x_tilde';

figure
subplot(4, 1, 1)
plot(t, x(:, 1), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 1), '--', 'LineWidth', 2)

subplot(4, 1, 2)
plot(t, x(:, 2), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 2), '--', 'LineWidth', 2)

subplot(4, 1, 3)
plot(t, x(:, 3), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 3), '--', 'LineWidth', 2)

subplot(4, 1, 4)
plot(t, x(:, 4), 'LineWidth', 2)
hold on
plot(t, x_bar(:, 4), '--', 'LineWidth', 2)

figure
plot(t, theta, 'LineWidth', 2)
hold on
plot(t, theta_bar, '--', 'LineWidth', 2)
% figure
% plot(t, x_tilde)



