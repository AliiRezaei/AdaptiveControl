function bmu_index = find_bmu(input_vector, weights)
    % Find the Best Matching Unit (BMU)
    min_dist = inf;
    bmu_index = 0;
    for j = 1:size(weights, 1)
        weight_vector = weights(j, :);
        dist = norm(input_vector - weight_vector);
        if dist < min_dist
            min_dist = dist;
            bmu_index = j;
        end
    end
end
