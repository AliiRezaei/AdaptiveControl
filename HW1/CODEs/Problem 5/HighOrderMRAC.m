function dx = HighOrderMRAC(t, states, r, problem)
    % Model Reference Adaptive Control of linear High Order systems (Direct method case).
    % In this case we have a linear sys as follows :
    %   xp' = Ap * xp + Bp * u
    % Where : 
    %   xp --> plant states        (n by 1)
    %   Ap --> plant system matrix (n by n)
    %   Bp --> plant input matrix  (n by 1)
    %   u  --> control signal      (1 by 1)
    % We want to design a control signal u such that :
    %    "the plant tracks the reference model"
    % Reference Model declare as below :
    %   xm' = Am * xm + Bm * r  (reference model)
    % Where : 
    %   xm --> reference model states         (n by 1)
    %   Am --> reference model system matrix  (n by n)
    %   Bm --> reference model input matrix   (n by 1)
    %   r  --> reference signal               (1 by 1)
    % Approache :
    %   In direct case, we choose the control signal u, such as :
    %       u = theta * xp + q * r
    %       then, we should calculate the best theta and q.
    % Inputs :
    %   t --> time
    %   states --> agumented states
    %   r --> reference signal
    %   problem --> contains plant and reference model params (problem formulation)
    %       problem.plant.Ap ~ Ap
    %       problem.plant.Bp ~ Bp
    %       problem.refModel.Am ~ Am
    %       problem.refModel.Bm ~ Bm
    %       problem.refModel.P  ~ P  (lyapanuve P matrix)
    %       problem.refModel.gamma ~ gamma (adaptation rate)
    %   Note : Am, P, Q satisfy following lyapunov equation :
    %       << Am^T * P + P * Am = - Q >> (Q declared in main program)
    % Outputs :
    %   dx --> derivative of agumented states
    
    % Actual system parameters :
    Ap = problem.plant.Ap;
    Bp = problem.plant.Bp;
    n = size(Ap, 1);
    nStates = numel(states);

    % Reference model parameters :
    Am = problem.refModel.Am;
    Bm = problem.refModel.Bm;
    P  = problem.refModel.P;
    gamma = problem.refModel.gamma;

    % Unpack states :
    xp_idx = 1:n; % measured indices
    xp = states(xp_idx); % plant measured states

    xm_idx = n+1:2*n; % ref model state indices
    xm  = states(xm_idx); % ref model states

    theta_hat_idx = 2*n+1:3*n; % controller param indices
    theta_hat = states(theta_hat_idx); % controller param

    % Designing control law :
    u = theta_hat' * xp + r(t);
    et = xp - xm; % tracking error

    % Derivative of states :
    dx = zeros(nStates, 1);
    dx(xp_idx)   = Ap * xp + Bp * u;
    dx(xm_idx)   = Am * xm + Bm * r(t);
    dx(theta_hat_idx) = - gamma * Bm' * P * et * xp';
    
end