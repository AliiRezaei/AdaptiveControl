function dx = MRAC(t, states, r, problem)

    % Note :
    %   MRAC implementation for any system with n* > 2 and known kp.
   
    % Plant parameters :
    Ap = problem.plant.Ap;
    Bp = problem.plant.Bp;
    Cp = problem.plant.Cp;
    kp = problem.plant.kp;
    n = size(Ap, 1);
    nStates = numel(states);

    % Refrence model parameters :
    Am = problem.refModel.Am;
    Bm = problem.refModel.Bm;
    Cm = problem.refModel.Cm;
    km = problem.refModel.km;
    m = size(Am, 1);

    % Designing parameters :
    LAMBDA = problem.designParam.LAMBDA;
    l = problem.designParam.l;
    gamma = problem.designParam.gamma;

    % Unpack states :
    xp_idx = 1:n;        % plant state indices
    xp = states(xp_idx); % plant states
    yp = Cp * xp;        % plant output
    
    xm_idx = n+1:n+m;    % ref model state indices
    xm = states(xm_idx); % ref model states
    ym = Cm * xm;        % ref model output
    
    w1_idx = n+m+1:2*n+m-1;
    w1 = states(w1_idx);
    
    w2_idx = 2*n+m:3*n+m-2;
    w2 = states(w2_idx);

    phi_bar_idx = 3*n+m-1:5*n+m-3;
    phi_bar = states(phi_bar_idx);

    e1_tilde_idx = 5*n+m-2:5*n+2*m-3;
    e1_tilde = states(e1_tilde_idx);

    zeta_tilde_idx = 5*n+2*m-2:nStates;
    zeta_tilde = states(zeta_tilde_idx);
    zeta_tilde = reshape(zeta_tilde, m, numel(zeta_tilde)/m);
    
    % Designing controller :
    w_bar = [w1; yp; w2];
    dzeta_tilde = [];
    zeta_bar = zeros(m, 1);
    for i = 1:numel(w_bar)
        tmp = Am * zeta_tilde(:, i) + Bm * w_bar(i);
        dzeta_tilde = [dzeta_tilde; tmp];
        zeta_bar(i) = Cm * tmp;
    end
    eps1 = phi_bar' * zeta_bar;
    u = phi_bar' * w_bar + kp / km * r(t);

    % Derivative of states :
    dx = zeros(nStates, 1);

    dx(xp_idx) = Ap * xp + Bp * u;
    dx(xm_idx) = Am * xm + Bm * r(t);

    dx(w1_idx) = LAMBDA * w1 + [zeros(n-2, 1); l] * u;
    dx(w2_idx) = LAMBDA * w2 + [zeros(n-2, 1); l] * yp;

    % Adaptation law :
    dx(phi_bar_idx) = - gamma * eps1 * zeta_bar;
    dx(zeta_tilde_idx) = dzeta_tilde;
    tmp = Am * e1_tilde + Bm * phi_bar' * w_bar;
    dx(e1_tilde_idx) = Cm * tmp;

end