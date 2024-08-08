function dx = NNAC(t, states, problem)
    % Adaptive Identification of linear High Order SISO systems.
    % In this case we have a linear sys as follows :
    %   xp' = Ap * xp + Bp * u
    % Where : 
    %   xp --> plant states        (n by 1)
    %   Ap --> plant system matrix (n by n)
    %   Bp --> plant input matrix  (n by 1)
    %   u  --> sys input           (1 by 1)
    % We want to estimate the system matrices (Ap and Bp), also system states.
    % This code use the Model (II) for parameterization.
    % Inputs :
    %   t --> time
    %   states --> agumented states
    %   u --> sys input
    %   problem --> contains plant and adaptive sys params (problem formulation)
    %       problem.plant.Ap ~ Ap
    %       problem.plant.Bp ~ Bp
    %       problem.adapt.Am ~ Am
    %       problem.adapt.P  ~ P  (lyapunov P matrix)
    %       problem.adapt.gamma ~ gamma
    %   Note : Am, P, Q satisfy following lyapunov equation :
    %       << Am^T * P + P * Am = - Q >> (Q declared in main program)
    % Outputs :
    %   dx --> derivative of agumented states
    
    % plant :
    f0 = problem.plant.fx; % f(x)
    g  = problem.plant.gx; % g(x)
    d  = problem.plant.d;  % d(t) disturbance
    
    % desired states :
    x_d = problem.desParam.x_d;

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
    x_idx = 1:n;                              % plant states indices
    x = states(x_idx);                        % plant states

    V_hat_idx = n+1:n+nInput*nHidden;         % input layer wights indices
    V_hat = states(V_hat_idx);                % input layer wights
    V_hat = reshape(V_hat, nHidden, nInput);  % reshape input layer wights

    W_hat_idx = n+nInput*nHidden+1:nStates;   % output layer wights indices
    W_hat = states(W_hat_idx);                % output layer wights
    W_hat = reshape(W_hat, nOutput, nHidden); % reshape output layer wights
    
    % init derivative of states :
    dx = zeros(nStates,1);

    % designing controller :
    e    = x - x_d(t)';                 % tracking error
    x_NN = [e; x_d(t)'];                % NN inputs
    u = - W_hat * kernel(V_hat * x_NN); % control signal

    % plant derivative :
    dx(x_idx)     = f0(x) + g(x) * u + d(t);

    % input layer weights adaption law :
    dV_hat = - eta(1) * (e' * inv(Ac) * ones(n, 1) * W_hat * (eye(nOutput) - diag((kernel(V_hat * x_NN)).^2)))' * x_NN' - ...
        gamma(1) * sqrt(e' * e) * V_hat; %#ok
    dx(V_hat_idx) = dV_hat(:);

    % output layer weights adaption law :
    dW_hat = - eta(2) * e' * inv(Ac) * ones(n, 1) * (kernel(V_hat * x_NN))' - ...
        gamma(2) * sqrt(e' * e) * W_hat; %#ok
    dx(W_hat_idx) = dW_hat;

end