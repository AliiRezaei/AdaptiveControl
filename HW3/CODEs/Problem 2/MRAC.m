function dx = MRAC(t, states, r, problem)
   
    % Plant parameters :
    Ap = problem.plant.Ap;
    Bp = problem.plant.Bp;
    Cp = problem.plant.Cp;
    n = size(Ap, 1);
    nStates = numel(states);

    % Refrence model parameters :
    Am = problem.refModel.Am;
    Bm = problem.refModel.Bm;
    Cm = problem.refModel.Cm;
    m = size(Am, 1);

    % Designing parameters :
    LAMBDA = problem.designParam.LAMBDA;
    l = problem.designParam.l;
    a = problem.designParam.a;
    gamma = problem.designParam.gamma;

    % Unpack states :
    xp_idx = 1:n; % plant state indices
    xp = states(xp_idx); % plant states
    yp = Cp * xp; % plant output
    
    xm_idx = n+1:n+m; % ref model state indices
    xm = states(xm_idx); % ref model states
    ym = Cm * xm; % ref model output
    
    w1_idx = n+m+1:2*n+m-1;
    w1 = states(w1_idx);
    
    w2_idx = 2*n+m:3*n+m-2;
    w2 = states(w2_idx);

    theta_hat_idx = 3*n+m-1:5*n+m-2;
    theta_hat = states(theta_hat_idx);

    w_bar_idx = 5*n+m-1:nStates;
    w_bar = states(w_bar_idx);
    
    % Designing controller :
    w = [r(t); w1; yp; w2];
    e1 = yp - ym;
    sign_kp = 1;
    dtheta_hat = - gamma * sign_kp * e1 .* [r(t); w_bar(2:n); yp; w_bar(end-n+2:end)];
    u = theta_hat' * w + dtheta_hat' * w_bar;

    % Derivative of states :
    dx = zeros(nStates,1);

    dx(xp_idx) = Ap * xp + Bp * u;
    dx(xm_idx) = Am * xm + Bm * r(t);

    dx(w1_idx) = LAMBDA * w1 + [zeros(n-2, 1); l] * u;
    dx(w2_idx) = LAMBDA * w2 + [zeros(n-2, 1); l] * yp;

    % Adaptation law :
    dx(theta_hat_idx) = dtheta_hat;

    % Low pass filter :
    dx(w_bar_idx) = - a * w_bar + w;


end