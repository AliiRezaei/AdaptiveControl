function dx = SOM(t, states, problem)
    
    % plant :
    f = problem.plant.fx;  % f(x)
    g  = problem.plant.gx; % g(x)
    d  = problem.plant.d;  % d(t) disturbance
    
    % desired states :
    x_d = problem.desParam.x_d;

    % NN dimensions :
    input_dim  = problem.SOM.InputDimension;
    input_neurons = prod(input_dim);

    output_dim  = problem.SOM.OutputDimension;
    output_neurons = prod(output_dim);
    
    nStates = numel(states);               % n augmented states
    n       = 2;                           % system order

    % SOM initial params :
    initial_radius = problem.SOM.initial_radius;
    initial_learning_rate = problem.SOM.initial_learning_rate;
    max_t = problem.SOM.max_epoch;

    % update SOM params :
    learning_rate = initial_learning_rate * (1 - t / (max_t + 1));
    radius = initial_radius * (1 - t / (max_t + 1));
    % learning_rate = initial_learning_rate * (0.001 / initial_learning_rate) ^ (t / max_t);
    % radius = initial_radius * (0.1 / initial_radius) ^ (t / max_t);

    % unpack states :
    x_idx = 1:n;                              % plant states indices
    x = states(x_idx);                        % plant states

    input_weights_idx = n+1:n+input_neurons*2*n;   % output layer wights indices
    input_weights = states(input_weights_idx);                % output layer wights
    input_weights = reshape(input_weights, input_neurons, 2*n); % reshape output layer wights

    output_weights_idx = n+input_neurons*2*n+1:nStates;   % output layer wights indices
    output_weights = states(output_weights_idx);                % output layer wights
    output_weights = reshape(output_weights, output_neurons, 1); % reshape output layer wights

    % input_weights_idx = n+1:n+input_neurons*n;   % output layer wights indices
    % input_weights = states(input_weights_idx);                % output layer wights
    % input_weights = reshape(input_weights, input_neurons, n); % reshape output layer wights
    % 
    % output_weights_idx = n+input_neurons*n+1:nStates;   % output layer wights indices
    % output_weights = states(output_weights_idx);                % output layer wights
    % output_weights = reshape(output_weights, output_neurons, 1); % reshape output layer wights

    % find input layer bmu :
    % input_vector = x';
    e = x - x_d(t)';
    input_vector = [e', x_d(t)];
    % input_vector = x_d(t);
    input_bmu_index = find_bmu(input_vector, input_weights);

    % get the row and column indices of the BMU
    [input_bmu_row, input_bmu_col] = ind2sub(input_dim, input_bmu_index);
    input_bmu = [input_bmu_row, input_bmu_col];


    % find output layer bmu :
    % output_vector = norm(input_weights(input_bmu));
    % output_vector = exp(-norm(input_weights(input_bmu))^2 / (2 * initial_radius^2));
    K = [5 3];
    output_vector = - K * x;
    output_bmu_index = find_bmu(output_vector, output_weights);

    % Get the row and column indices of the BMU
    [output_bmu_row, output_bmu_col] = ind2sub(output_dim, output_bmu_index);
    output_bmu = [output_bmu_row, output_bmu_col];
    
    % update input layer weights :
    dinput_weights = update_weights(input_weights, input_vector, input_dim, input_bmu, radius, learning_rate);

    % update output layer weights :
    doutput_weights = update_weights(output_weights, output_vector, output_dim, output_bmu, radius, learning_rate);

    % init derivative of states :    
    dx = zeros(nStates,1);

    % plant derivative :
    % u = output_vector;
    u = output_weights(output_bmu_row, output_bmu_col);
    dx(x_idx)     = f(x) + g(x) * u + d(t);

    % input layer weights update rule :
    dx(input_weights_idx) = dinput_weights(:);

    % output layer weights update rule :
    dx(output_weights_idx) = doutput_weights(:);

end