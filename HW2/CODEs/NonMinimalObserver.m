function dx = NonMinimalObserver(t, states, u, problem)
   
    % Plant parameters :
    Ap = problem.plant.Ap;
    Bp = problem.plant.Bp;
    Cp = problem.plant.Cp;
    n = size(Ap, 1);
    nStates = numel(states);

    % Designing parameters :
    lambda = problem.designParam.lambda;
    LAMBDA = problem.designParam.LAMBDA;
    l = problem.designParam.l;
    gamma = problem.designParam.gamma;

    % Unpack states :
    x_idx = 1:n; % measured indices
    x = states(x_idx); % measured states
    y = Cp * x;

    x1_hat_idx = n+1;
    x1_hat = states(x1_hat_idx);
    y_hat = x1_hat;
    
    w1_hat_idx = n+2:2*n;
    w1_hat = states(w1_hat_idx);
    
    w2_hat_idx = 2*n+1:3*n-1;
    w2_hat = states(w2_hat_idx);

    theta_hat_idx = 3*n:nStates;
    theta_hat = states(theta_hat_idx);

    w_hat = [u(t); w1_hat; y; w2_hat];

    % Derivative of states :
    dx = zeros(nStates,1);
    dx(x_idx)     = Ap * x + Bp * u(t);

    % Non-Minimal realization :
    dx(x1_hat_idx) = - lambda * x1_hat + theta_hat' * w_hat;
    dx(w1_hat_idx) =   LAMBDA * w1_hat + [zeros(n-2, 1); l] * u(t);
    dx(w2_hat_idx) =   LAMBDA * w2_hat + [zeros(n-2, 1); l] * y;

    % Adaptation law :
    e = y_hat - y;
    dx(theta_hat_idx) = - gamma * e * w_hat;

end