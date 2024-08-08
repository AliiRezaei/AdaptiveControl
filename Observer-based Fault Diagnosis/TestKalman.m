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

y = zeros(4, nt);
y(:, 1) = x(:, 1);

x_bar = zeros(4, nt);
x_bar(:, 1) = randn(4, 1);

x_tilde = zeros(4, nt);
x_tilde(:, 1) = x(:, 1) - x_bar(:, 1);

theta = zeros(nt, 1);
theta_bar = zeros(nt, 1);
theta_tilde = zeros(nt, 1);

P = 2 * eye(4);
R = 10;
Q = 1 * eye(4);

u = ones(nt, 1);

for kk = 1:nt - 1

    % u(kk) = 1;
    % u(kk) = B \ ((A_tilde - A) * x(:, kk) - f(x(3, kk))) / (1 - theta(kk));
    PI = - B * u(kk);

    if kk*dt > 2
        theta(kk+1) = 1.5;
    end

    w = 0 * randn(4, 1);
    v = 0 * randn;
    x(:, kk+1) = A * x(:, kk) + B * u(kk) + w;
    y(:, kk)   = x(:, kk) + v;
    K = P * R^(-1);
    x_bar(:, kk+1) = (1 - K) *(A * x_bar(:, kk) + B * u(kk)) + K * y(:, kk);
    P = (eye(4) - K) * (A * P * A' + Q);
    if any(eig(P) < 0)
        kk
        eig(P)
        break;
    end

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
figure
plot(t, x_tilde)



