clc
clear
% close all
set(0, 'defaultTextInterpreter', 'latex');

%% Problem Definition

% plant :
problem.plant.fx = @(x) [x(2); -5*x(1)^3 - 2*x(2)]; % f(x)
problem.plant.gx = @(x) [0; 1];                           % g(x)
problem.plant.d  = @(t) 0 * (t > 10);                     % d(t) disturbance
n = 2;

% desired states :
problem.desParam.x_d = @(t) [sin(t), cos(t)];

% constants (for controller designing) :
problem.desParam.Ac    = diag([-20, -20]);
problem.desParam.eta   = [25; 25];
problem.desParam.gamma = [0.001; 0.001];

%% Neural Network Params

% layers :
inputLayerSize  = 2 * n;                  % number of neurons in input  layer
hiddenLayerSize = 2 * inputLayerSize + 1; % number of neurons in hidden layer
outputLayerSize = 1;                      % number of neurons in output layer
problem.NN.inputLayerSize  = inputLayerSize;
problem.NN.hiddenLayerSize = hiddenLayerSize;
problem.NN.outputLayerSize = outputLayerSize;

% % RBF kernel function :
% mu    = zeros(hiddenLayerSize, outputLayerSize); % mean
% sigma = ones(hiddenLayerSize,  outputLayerSize); % var
% problem.NN.kernel = @(phi) exp((phi - mu) ./ sigma.^2);

% tanh kernel function :
problem.NN.kernel = @(phi) tanh(phi);

% weights (initialization):
V = 0.1 * randn(hiddenLayerSize, inputLayerSize);  % hidden layer weights
W = 0.1 * randn(outputLayerSize, hiddenLayerSize); % output weights

%% Simulate System

dt = 0.05;               % time step [seconds]
SimTime = 50;            % maximum simulation time [seconds]
t = (0:dt:SimTime)';     % time span
x_0  = [0.5; 0.5];       % plant init cond
V_hat_0 = V;             % hidden layer weights init cond
W_hat_0 = W;             % output weights init cond

InitCond = [x_0; V_hat_0(:); W_hat_0(:)]; % initial condition
odeFunc  = @(t, x) NNAC(t, x, problem);   % ode function
[~, X]   = ode45(odeFunc, t, InitCond);   % solve ode
nStates  = size(X, 1);                    % number of states

%% Unpack States

% desired states :
x_d = problem.desParam.x_d(t);

% plant states :
x   = X(:, 1:n);

% input layer weights :
c = 1;
V_hat = zeros(hiddenLayerSize, inputLayerSize, nStates);
for i = 1:hiddenLayerSize
    for j = 1:inputLayerSize
        V_hat(i, j, :) = X(:, n+c);
        c = c + 1;
    end
end

% output layer weights :
W_hat = zeros(outputLayerSize, hiddenLayerSize, nStates);
for i = 1:outputLayerSize
    for j = 1:hiddenLayerSize
        W_hat(i, j, :) = X(:, n+c);
        c = c + 1;
    end
end

% control signal :
u = zeros(nStates, 1);
for k = 1:nStates
    u(k, 1) = - W_hat(:, :, k) * problem.NN.kernel(V_hat(:, :, k) * [(x(k, :) - x_d(k, :))'; x_d(k, :)']);
end

%% Trained Network

net = @(x) - W_hat(:, :, end) * problem.NN.kernel(V_hat(:, :, end) * x);

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
for i = 1:hiddenLayerSize
    for j = 1:inputLayerSize
        plot(t, squeeze(V_hat(i, j, :)), 'LineWidth', 2)
        hold on
    end
end
ylabel('$\hat{V}(t)$', 'FontSize', 15)
title('Input Layer Weights', 'FontSize', 15)

subplot(2,1,2)
for i = 1:outputLayerSize
    for j = 1:hiddenLayerSize
        plot(t, squeeze(W_hat(i, j, :)), 'LineWidth', 2)
        hold on
    end
end
xlabel('$t$ $[sec]$', 'FontSize', 15)
ylabel('$\hat{W}(t)$', 'FontSize', 15)
title('Output Layer Weights', 'FontSize', 15)

% plot control signal :
figure
plot(t, u, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('u(t)', 'FontSize', 15)
title('Control Signal', 'FontSize', 15)

% plot tracking error :
figure
plot(t, x - x_d, 'LineWidth', 2)
xlabel('t [sec]', 'FontSize', 15)
ylabel('e(t)', 'FontSize', 15)
title('Tracking Error', 'FontSize', 15)
