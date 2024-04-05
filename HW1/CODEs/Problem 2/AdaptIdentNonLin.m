function dx = AdaptIdentNonLin(t, states, u, problem)
    % Adaptive Identification of Nonlinear scalar systems.
    % In this case we have a nonlinear sys as follows :
    %   xp' = ap * xp + alpha * f(xp) + kp * g(u)
    % Note that f(xp) and g(u) are known functions.
    % Where : 
    %   xp --> sys state               (1 by 1)
    %   ap --> linear term gain        (1 by 1)
    %   alpha --> nonlinear term gain  (1 by 1)
    %   f(xp) --> nonlinear scalar function of xp
    %   kp --> feedforward gain        (1 by 1)
    %   u  --> sys input               (1 by 1)
    %   g(u) --> nonlinear scalar function of u
    % We want to estimate the system params (ap, alpha, kp), also system state.
    % This code cover the both Model (I) and Model (II) parameterization.
    % Model (I) :
    %   xp_hat' = ap_hat * xp_hat + alpha_hat * f(xp) + kp_hat * g(u)
    % Model (II) :
    %   xp_hat' = am * xp_hat + (ap_hat - am) * xp + alpha_hat * f(xp) + kp_hat * g(u)
    % For using Model (I)  --> set the problem.adapt.modelID to 1
    % For using Model (II) --> set the problem.adapt.modelID to 2
    % Inputs :
    %   t --> time
    %   states --> agumented states
    %   u --> sys input
    %   problem --> contains plant and adaptive sys params (problem formulation)
    %       problem.plant.ap ~ ap
    %       problem.plant.kp ~ kp
    %       problem.plant.alpha ~ alpha
    %       problem.adapt.am ~ am
    %       problem.adapt.gamma ~ gamma
    %       problem.adapt.f ~ f(xp)
    %       problem.adapt.g ~ g(u)
    %       problem.adapt.modelID ~ modelID
    % Outputs :
    %   dx --> derivative of agumented states
    
    % Unpack states :
    xp        = states(1); % system measured state          derivative --> dx(1)
    xp_hat    = states(2); % system estimated state         derivative --> dx(2)
    ap_hat    = states(3); % estimated linear term gain     derivative --> dx(3)
    alpha_hat = states(4); % estimated nonlinear term gain  derivative --> dx(4)
    kp_hat    = states(5); % estimated feedforward gain     derivative --> dx(5)

    % Real system parameters :
    %   Note : xp' = ap * xp + alpha * f(xp) + kp * g(u)
    ap    = problem.plant.ap;     % linear term gain
    alpha = problem.plant.alpha;  % nonlinear term gain
    kp    = problem.plant.kp;     % actual feedforward gain
    f     = problem.plant.f;      % nonlinear function of xp [f(xp)]
    g     = problem.plant.g;      % nonlinear function of  u [g(u)]

    % Adaptation law parameters :
    am      =  problem.adapt.am;     % hurwitz param
    gamma   =  problem.adapt.gamma;  % adaptation rate
    modelID = problem.adapt.modelID; % model ID
    e = xp_hat - xp;              % estimation error

    % Derivative of states (using Model (I) case):
    if modelID == 1
        dx    = zeros(5, 1);
        dx(1) = ap * xp + alpha * f(xp) + kp * g(u(t));
        dx(2) = ap_hat * xp_hat + alpha_hat * f(xp) + kp_hat * g(u(t));
        dx(3) = - gamma * e * xp;
        dx(4) = - gamma * e * f(xp);
        dx(5) = - gamma * e * g(u(t));
    end

    % Derivative of states (using Model (II) case):
    if modelID == 2
        dx    = zeros(5, 1);
        dx(1) = ap * xp + alpha * f(xp) + kp * g(u(t));
        dx(2) = am * xp_hat + (ap_hat - am) * xp + alpha_hat * f(xp) + kp_hat * g(u(t));
        dx(3) = - gamma * e * xp;
        dx(4) = - gamma * e * f(xp);
        dx(5) = - gamma * e * g(u(t));
    end
        
end