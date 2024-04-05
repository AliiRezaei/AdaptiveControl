function dx = AdaptIdent(t, states, u, modelID)
    % Adaptive Identification of linear scaler systems.
    % In this case we have a linear sys as follows :
    %   xp' = ap * xp + kp * u
    % Where : 
    %   xp --> sys state        (1 by 1)
    %   ap --> feedback gain    (1 by 1)
    %   kp --> feedforward gain (1 by 1)
    %   u  --> sys input        (1 by 1)
    % We want to estimate the system params (ap and kp), also system state.
    % This code cover the both Model (I) and Model (II) parameterization.
    % Model (I) :
    %   xp_hat' = ap_hat * xp_hat + kp_hat * u
    % Model (II) :
    %   xp_hat' = ap_hat * xp_hat + kp_hat * u
    % For using Model (I)  --> set the modelID to 1
    % For using Model (II) --> set the modelID to 2
    % Inputs :
    %   t --> time
    %   states --> agumented states
    %   u --> sys input
    %   modelID --> 1 for using Model (I) or 2 for using Model (II)
    % Outputs :
    %   dx --> derivative of agumented states
    
    % Unpack states :
    xp     = states(1); % system measured state       derivative --> dx(1)
    xp_hat = states(2); % system estimated state      derivative --> dx(2)
    ap_hat = states(3); % estimated feedback gain     derivative --> dx(3)
    kp_hat = states(4); % estimated feedforward gain  derivative --> dx(4)

    % Real system parameters :
    %   Note : xp' = ap * xp + kp * u
    ap = -1; % actual feedback gain
    kp =  1; % actual feedforward gain

    % Adaptation law parameters :
    am    =  -2;         % hurwitz param
    gamma =  10;         % adaptation rate
    e = xp_hat - xp;     % estimation error

    % Derivative of states (using Model (I) case):
    if modelID == 1
        dx    = zeros(4, 1);
        dx(1) = ap * xp + kp * u(t);
        dx(2) = ap_hat * xp_hat + kp_hat * u(t);
        dx(3) = - gamma * e * xp_hat;
        dx(4) = - gamma * e * u(t);
    end

    % Derivative of states (using Model (II) case):
    if modelID == 2
        dx    = zeros(4, 1);
        dx(1) = ap * xp + kp * u(t);
        dx(2) = am * xp_hat + (ap_hat - am) * xp + kp_hat * u(t);
        dx(3) = - gamma * e * xp;
        dx(4) = - gamma * e * u(t);
    end

end