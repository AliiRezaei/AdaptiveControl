function dweights = update_weights(weights, x, grid, bmu, radius, learning_rate)
    dweights = zeros(size(weights));
    % Update the weights of the BMU and its neighbors
    for j = 1:size(weights, 1)
        % Get the row and column indices of the current neuron
        [neuron_row, neuron_col] = ind2sub(grid, j);
        distance_to_bmu = norm([neuron_row, neuron_col] - bmu);
        if distance_to_bmu <= radius
            % Calculate the influence based on the neighborhood function
            influence = exp(-distance_to_bmu^2 / (2 * radius^2));
            % Update weights
            dweights(j, :) = learning_rate * influence * (x - weights(j, :));
        end
    end
end