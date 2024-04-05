function dx = DirectMRAC(t, states, r, problem)
    % Model Reference Adaptive Control Direct method.
    % In this case we have a linear sys as follows :
    %   xp' = ap * xp + kp * u  (plant)
    % Where : 
    %   xp --> plant state            (1 by 1)
    %   ap --> plant feedback gain    (1 by 1)
    %   kp --> plant feedforward gain (1 by 1)
    %   u  --> control signal         (1 by 1)
    % We want to design a control signal u such that :
    %    "the plant tracks the reference model"
    % Reference Model declare as below :
    %   xm' = am * xm + km * r  (reference model)
    % Where : 
    %   xm --> reference model state            (1 by 1)
    %   am --> reference model feedback gain    (1 by 1)
    %   km --> reference model feedforward gain (1 by 1)
    %   r  --> reference signal                 (1 by 1)
    % Approache :
    %   In direct case, we choose the control signal u, such as :
    %       u = theta * xp + k * r
    %       then, we should calculate the best theta and k.
    % Inputs :
    %   t --> time
    %   states --> agumented states
    %   r --> reference signal
    % Outputs :
    %   dx --> derivative of agumented states
    
    % Unpack states :
    xp        = states(1); % plant measured state  derivative --> dx(1)
    xm        = states(2); % ref model state       derivative --> dx(2)
    theta_hat = states(3); % controller param      derivative --> dx(3)
    k_hat     = states(4); % controller param      derivative --> dx(4)

    % Real system parameters (plant) :
    %   Note : xp' = ap * xp + kp * u
    ap = problem.plant.ap; % actual sys feedback gain
    kp = problem.plant.kp; % actual sys feedforward gain

    % Reference system parameters :
    %   Note : xm' = am * xm + km * r
    am = problem.refModel.am; % ref model feedback gain
    km = problem.refModel.km; % ref model feedforward gain
    gamma = problem.refModel.gamma; % adaptation rate

    % Designing control law :
    u = theta_hat * xp + k_hat * r(t);
    et = xp - xm; % tracking error

    % Derivative of states (direct method case):
    dx    = zeros(4, 1);
    dx(1) = ap * xp + kp * u;
    dx(2) = am * xm + km * r(t);
    dx(3) = - gamma * sign(kp) * et * xp;
    dx(4) = - gamma * sign(kp) * et * r(t);

end