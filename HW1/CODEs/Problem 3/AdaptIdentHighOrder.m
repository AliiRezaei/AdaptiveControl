function dx = AdaptIdentHighOrder(t, states, u, problem)
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
    %   actulaSys --> actula system (plant) params. 
    %   contains :
    %       actulaSys.A ~ Ap
    %       actulaSys.B ~ Bp
    %   adaptSys --> adaptive system params. 
    %   contains :
    %       adaptSys.Am ~ Hurwitz matrix with size Ap
    %       adaptSys.P  ~ lyapunov P matrix
    %       adaptSys.Q  ~ lyapunov Q matrix
    %       adaptSys.gamma  ~ adaptation rate
    %   Note : Am, P, Q satisfy following lyapunov equation :
    %       << Am^T * P + P * Am = -Q >>
    % Outputs :
    %   dx --> derivative of agumented states
    
    % actual system parameters :
    Ap = problem.plant.Ap;
    Bp = problem.plant.Bp;
    n = size(Ap, 1);
    nStates = numel(states);

    % adaptive system parameters :
    Am = problem.adapt.Am;
    P  = problem.adapt.P;
    gamma = problem.adapt.gamma;

    % Unpack states :
    x_idx = 1:n; % measured indices
    x = states(x_idx); % measured states

    x_hat_idx = n+1:2*n; % estimation indices
    x_hat  = states(x_hat_idx); % estimated states

    A_hat_idx = 2*n+1:2*n+n^2;  % system matrix indices
    A_hat = states(A_hat_idx);  % estimated system matrix
    A_hat = reshape(A_hat,n,n); % reshape estimated sys matrix (n by n) 

    B_hat_idx = 2*n+1+n^2:nStates;  % system input matrix indices
    B_hat = states(B_hat_idx);      % estimated input matrix
    B_hat = reshape(B_hat,n,1);     % reshape estimated input matrix (n by 1)

    % Derivative of states :
    e  = x_hat - x;
    dx = zeros(nStates,1);
    dx(x_idx)     = Ap * x + Bp * u(t);
    dx(x_hat_idx) = Am * x_hat + (A_hat - Am) * x + B_hat * u(t);
    dx(A_hat_idx) = - gamma * P * e * x';
    dx(B_hat_idx) = - gamma * P * e * u(t)';

end