clc
clear
close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% plant :
problem.plant.fx = @(x) [x(2); - x(2) - x(1)]; % f(x)
problem.plant.gx = @(x) [0; 1];                           % g(x)
problem.plant.d  = @(t) 0 * (t > 10);                     % d(t) disturbance
n = 2;

% desired states :
% problem.desParam.x_d = @(t) 0.1 * [sin(10 * t), cos(10 * t)];
% problem.desParam.x_d = @(t) [sin(0.5*t) + cos(0.25*t), -0.25*sin(0.25*t) + 0.5*cos(0.5*t)];
% problem.desParam.x_d = @(t) 0 .* [sin(0.5*t) + cos(0.25*t), -0.25*sin(0.25*t) + 0.5*cos(0.5*t)];
problem.desParam.x_d = @(t) .1 .* [1, 0];

%% Neural Network Params

% layers :
input_grid_rows = 10; % Number of rows in SOM grid
input_grid_cols = 1; % Number of columns in SOM grid
problem.SOM.InputDimension  = [input_grid_rows, input_grid_cols];

output_grid_rows = input_grid_rows * input_grid_cols; % Number of rows in SOM grid
output_grid_cols = 1; % Number of columns in SOM grid
problem.SOM.OutputDimension  = [output_grid_rows, output_grid_cols];

input_neurons  = input_grid_rows  * input_grid_cols;
output_neurons = output_grid_rows * output_grid_cols;

% weights (initialization):
input_weights  = 0.1 * randn(input_neurons,  2*n); % input  weights
output_weights = 0.1 * randn(output_neurons, 1); % output weights

%% Simulate System

dt = 0.05;               % time step [seconds]
SimTime = 100;            % maximum simulation time [seconds]
t = (0:dt:SimTime)';     % time span

problem.SOM.initial_learning_rate = 0.1; % Initial learning rate
% problem.SOM.initial_radius = 10; % Initial radius
problem.SOM.initial_radius = max(input_grid_rows, input_grid_cols) / 2; % Initial radius
problem.SOM.max_epoch = SimTime / dt; % Initial learning rate

x_0  = [0.5; 0.5];       % plant init cond

InitCond = [x_0; input_weights(:); output_weights(:)]; % initial condition
odeFunc  = @(t, x) SOM(t, x, problem);   % ode function
[~, X]   = ode45(odeFunc, t, InitCond);   % solve ode
nStates  = size(X, 1);                    % number of states

%% Unpack States

% desired states :
x_d = problem.desParam.x_d(t);

% plant states :
x   = X(:, 1:n);

% input layer weights :
c = 1;
input_weights = zeros(input_neurons, n, nStates);
for i = 1:input_neurons
    for j = 1:2*n
        input_weights(i, j, :) = X(:, n+c);
        c = c + 1;
    end
end

output_weights = zeros(output_neurons, 1, nStates);
for i = 1:output_neurons
    output_weights(i, 1, :) = X(:, n+c);
    c = c + 1;
end

% % control signal :
% u = zeros(nStates, 1);
% for k = 1:nStates
%     u(k, 1) = - W_hat(:, :, k) * problem.NN.kernel(V_hat(:, :, k) * [(x(k, :) - x_d(k, :))'; x_d(k, :)']);
% end

% %% Trained Network
% 
% net = @(x) - W_hat(:, :, end) * problem.NN.kernel(V_hat(:, :, end) * x);

%% Plots and Results

% plot x and x_d :
figure
subplot(2,1,1)
plot(t, x(:, 1), 'LineWidth', 2)
hold on
plot(t, x_d(:, 1), '--', 'LineWidth', 2)
ylabel('$x_1(t)$', 'FontSize', 15)
title('Plant Actual and Desired States', 'FontSize', 15)
legend('$x_1(t)$','$x_{1d}(t)$', 'interpreter', 'latex')

subplot(2, 1, 2)
plot(t, x(:, 2), 'LineWidth', 2)
hold on
plot(t, x_d(:, 2), '--', 'LineWidth', 2)
ylabel('$x_2(t)$', 'FontSize', 15)
legend('$x_2(t)$','$x_{2d}(t)$', 'interpreter', 'latex')
xlabel('$t$ $[sec]$', 'FontSize', 15)

% plot estimated NN weights :
figure
subplot(2,1,1)
Color = hsv(input_neurons);
for i = 1:input_neurons
    for k = 1:nStates
        center_x = squeeze(input_weights(i, 1, :));
        center_y = squeeze(input_weights(i, 2, :));
        plot(center_x, center_y, '-o', 'MarkerEdgeColor', Color(i, :), 'MarkerSize', 8, 'MarkerFaceColor', Color(i, :))
        hold on
    end
    Color(i, :) = 0.99 * Color(i, :);
end
ylabel('$\hat{V}(t)$', 'FontSize', 15)
title('Input Layer Weights', 'FontSize', 15)

subplot(2,1,2)
Color = hsv(output_neurons);
for i = 1:output_neurons
    for k = 1:nStates
        center = squeeze(output_weights(i, 1, :));
        plot(center, '-o', 'MarkerEdgeColor', Color(i, :), 'MarkerSize', 8, 'MarkerFaceColor', Color(i, :))
        hold on
    end
    Color(i, :) = 0.99 * Color(i, :);
end
xlabel('$t$ $[sec]$', 'FontSize', 15)
ylabel('$\hat{W}(t)$', 'FontSize', 15)
title('Output Layer Weights', 'FontSize', 15)

% % plot control signal :
% figure
% plot(t, u, 'LineWidth', 2)
% xlabel('t [sec]', 'FontSize', 15)
% ylabel('u(t)', 'FontSize', 15)
% title('Control Signal', 'FontSize', 15)

% plot tracking error :
figure
plot(t, x - x_d, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('e(t)', 'FontSize', 15)
title('Tracking Error', 'FontSize', 15)
