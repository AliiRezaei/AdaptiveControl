function dx = NNACIS(t, states, problem)
    % Neural Network-Based State Feedback Adaptive Control pf Affine Nonlinear Systems with Input Saturation.
    % This function implements a controller for a geometric nonlinear system defined as:
    %   x' = f(x) + g(x) * sat(u) + d(t)
    % Where: 
    %   x    --> Plant states                        (n by 1)
    %   f(x) --> Nonlinear smooth function of states (n by 1)
    %   g(x) --> Nonlinear smooth function of states (n by m)
    %   u    --> System input                        (m by 1)
    %   d    --> Disturbance                         (1 by 1)
    % "n" is the system order, and "m" is the number of system inputs.
    %
    % The goal is to design a Neural Network-Based Adaptive Controller
    % that drives the system states to the desired states without prior
    % knowledge of the system dynamics.
    %
    % This is achieved using a two-layer feedforward Neural Network
    % with "n" neurons in the input layer, "2n+1" neurons in the hidden layer,
    % and "m" neurons in the output layer.
    %
    % The controller is expressed as:
    %   u = -W * PHI(V * x_NN)
    % Where: 
    %   V    --> Input layer weights 
    %   W    --> Output layer weights 
    %   PHI  --> Network kernel function 
    %   x_NN --> NN inputs 
    %
    % The adaptation laws for V and W are derived accordingly.
    %
    % Inputs:
    %   t       --> Time
    %   states  --> Augmented states vector
    %   problem --> Structure containing plant, NN, and design parameters
    %       problem.plant.fx           --> Function handle for f(x)
    %       problem.plant.gx           --> Function handle for g(x)
    %       problem.plant.d            --> Function handle for disturbance d(t)
    %       problem.NN.inputLayerSize  --> Input layer size
    %       problem.NN.hiddenLayerSize --> Hidden layer size
    %       problem.NN.outputLayerSize --> Output layer size
    %       problem.desParam.Ac        --> Hurwitz matrix (n by n)
    %       problem.desParam.eta       --> Learning rate vector (2 by 1)
    %       problem.desParam.gamma     --> E-modification constants (2 by 1)
    %
    % Outputs:
    %   dx --> Derivative of augmented states
    
    % plant :
    f = problem.plant.fx;  % f(x)
    g  = problem.plant.gx; % g(x)
    d  = problem.plant.d;  % d(t) disturbance
    
    % desired states :
    x_d = problem.desParam.x_d;

    % auxiliary dynamics params :
    K = problem.auxiliary.K;
    h = problem.auxiliary.h;
    
    % input saturation bounds :
    u_min = problem.sat.u_min;
    u_max = problem.sat.u_max;

    % designing params :
    Ac    = problem.desParam.Ac;
    eta   = problem.desParam.eta;
    gamma = problem.desParam.gamma;

    % NN kernel function :
    kernel = problem.NN.kernel;

    % NN dimensions :
    nInput  = problem.NN.inputLayerSize;  % input  layer size
    nHidden = problem.NN.hiddenLayerSize; % hidden layer size
    nOutput = problem.NN.outputLayerSize; % output layer size
    nStates = numel(states);              % n augmented states
    n       = size(Ac, 1);                % system order

    % unpack states :
    x_idx = 1:n;                                % plant states indices
    x = states(x_idx);                          % plant states

    % unpack states :
    zeta_idx = n+1:2*n;                         % auxiliary states indices
    zeta = states(zeta_idx);                    % auxiliary states

    V_hat_idx = 2*n+1:2*n+nInput*nHidden;       % input layer wights indices
    V_hat = states(V_hat_idx);                  % input layer wights
    V_hat = reshape(V_hat, nHidden, nInput);    % reshape input layer wights

    W_hat_idx = 2*n+nInput*nHidden+1:nStates;   % output layer wights indices
    W_hat = states(W_hat_idx);                  % output layer wights
    W_hat = reshape(W_hat, nOutput, nHidden);   % reshape output layer wights
    
    % init derivative of states :
    dx = zeros(nStates,1);

    % designing controller :
    e       = x - x_d(t)';                    % tracking error
    e_tilde = e - zeta;                       % auxiliary tracking error
    x_NN    = [e_tilde; x_d(t)'];             % NN inputs
    u       = - W_hat * kernel(V_hat * x_NN); % control signal
    
    % saturated control signal :
    if u > u_max
        sat_u = u_max;
    elseif u < u_min
        sat_u = u_min;
    else
        sat_u = u;
    end

    % ctrl signal part that is not applied to the system due to actuator saturation :
    delta_u = u - sat_u;

    % plant derivative :
    dx(x_idx)     = f(x) + g(x) * sat_u + d(t);

    % auxiliary satates derivative :
    dx(zeta_idx)  = K * zeta + h * delta_u;

    % input layer weights adaption law :
    dV_hat = - eta(1) * (e_tilde' * (-inv(Ac) * ones(n, 1) + inv(K) * h * delta_u) * (-W_hat * (eye(nHidden ...
        ) - diag((kernel(V_hat * x_NN)).^2))))' * x_NN' - ...
        gamma(1) * sqrt(e' * e) * V_hat - gamma(1) * sqrt(zeta' * zeta) * V_hat; %#ok
    dx(V_hat_idx) = dV_hat(:);

    % output layer weights adaption law :
    dW_hat = - eta(2) * e_tilde' * (-inv(Ac) * ones(n, 1) + inv(K) * h * delta_u) * (-kernel(V_hat * x_NN))' - ...
        gamma(2) * sqrt(e' * e) * W_hat - gamma(2) * sqrt(zeta' * zeta) * W_hat; %#ok
    dx(W_hat_idx) = dW_hat(:);

end