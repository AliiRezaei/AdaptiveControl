function dx = IndirectMRAC(t, states, r, problem)
    % Model Reference Adaptive Control Indirect method.
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
    %   In indirect case, we choose the control signal u, such as :
    %       u = (am - ap) / (kp) * xp + (km) / (kp) * r
    %       but, ap and kp are unknown!
    %       so, we use ap_hat and kp_hat instead of ap and kp.
    % Inputs :
    %   t --> time
    %   states --> agumented states
    %   r --> reference signal
    % Outputs :
    %   dx --> derivative of agumented states
    
    % Unpack states :
    xp     = states(1); % plant measured state             derivative --> dx(1)
    xp_hat = states(2); % plant estimated state            derivative --> dx(2)
    xm     = states(3); % ref model state                  derivative --> dx(3)
    ap_hat = states(4); % plant estimated feedback gain    derivative --> dx(4)
    kp_hat = states(5); % plant estimated feedforward gain derivative --> dx(5)

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
    u = (am - ap_hat) / (kp_hat) * xp + (km) / (kp_hat) * r(t); % control law
    ei = xp_hat - xp; % identification error

    % Derivative of states (indirect method case):
    dx    = zeros(4, 1);
    dx(1) = ap * xp + kp * u;
    dx(2) = am * xp_hat + (ap_hat - am) * xp + kp_hat * u;
    dx(3) = am * xm + km * r(t);
    dx(4) = - gamma * ei * xp;
    dx(5) = - gamma * ei * u;

end