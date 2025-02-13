clc
clear
close all

% Load the simpleclass_dataset
% [x, t] = simpleclass_dataset;
[x, t] = iris_dataset;
data = x'; % Transpose for easier processing

% SOM Parameters
grid_rows = 4; % Number of rows in SOM grid
grid_cols = 1; % Number of columns in SOM grid
input_dim = size(data, 2); % Dimension of input data
max_iter = 100; % Maximum number of iterations
initial_learning_rate = 0.1; % Initial learning rate
initial_radius = max(grid_rows, grid_cols) / 2; % Initial radius

% Initialize the weights randomly in a 2D array
weights = rand(grid_rows * grid_cols, input_dim);

% Function to calculate Euclidean distance
euclidean_distance = @(vec1, vec2) sqrt(sum((vec1 - vec2).^2));

% Training the SOM
for iter = 1:max_iter
    % Decrease learning rate and radius
    learning_rate = initial_learning_rate * (1 - iter / (max_iter + 1));
    radius = initial_radius * (1 - iter / (max_iter + 1));
    
    for i = 1:size(data, 1)
        input_vector = data(i, :);
        
        % Find the Best Matching Unit (BMU)
        min_dist = inf;
        bmu_index = 0;
        for j = 1:size(weights, 1)
            weight_vector = weights(j, :);
            dist = euclidean_distance(input_vector, weight_vector);
            if dist < min_dist
                min_dist = dist;
                bmu_index = j;
            end
        end
        
        % Get the row and column indices of the BMU
        [bmu_row, bmu_col] = ind2sub([grid_rows, grid_cols], bmu_index);
        
        % Update the weights of the BMU and its neighbors
        for j = 1:size(weights, 1)
            % Get the row and column indices of the current neuron
            [neuron_row, neuron_col] = ind2sub([grid_rows, grid_cols], j);
            distance_to_bmu = norm([neuron_row, neuron_col] - [bmu_row, bmu_col]);
            if distance_to_bmu <= radius
                % Calculate the influence based on the neighborhood function
                influence = exp(-distance_to_bmu^2 / (2 * radius^2));
                % Update weights
                weights(j, :) = weights(j, :) + ...
                    learning_rate * influence * (input_vector - weights(j, :));
            end
        end
    end
end

% Plot the results
figure;
hold on;
scatter(data(:, 1), data(:, 2), 'filled');
for j = 1:size(weights, 1)
    % Get the row and column indices of the neuron for plotting
    [neuron_row, neuron_col] = ind2sub([grid_rows, grid_cols], j);
    plot(weights(j, 1), weights(j, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end
title('SOM Clustering of simpleclass\_dataset');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend('Data Points', 'Neuron Positions');
hold off;