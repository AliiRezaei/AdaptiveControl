function dx = MRAC2(t, states, r, problem)

    % Note :
    %   MRAC implementation for any system with n* > 2 and unknown kp.
   
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

    phi_idx = 3*n+m-1:5*n+m-3;
    phi = states(phi_idx);

    psi_idx = 5*n+m-2;
    psi = states(psi_idx);

    e2_virtual_idx = 5*n+m-1:5*n+2*m-2;
    e2_virtual = states(e2_virtual_idx);

    zeta_virtual_idx = 5*n+2*m-1:nStates;
    zeta_virtual = states(zeta_virtual_idx);
    zeta_virtual = reshape(zeta_virtual, m, numel(zeta_virtual)/m);
    
    % Designing controller :
    w_bar = [w1; yp; w2];
    zeta = (Cm * zeta_virtual)';
    % zeta = zeta / (1+zeta'*zeta);
    % e1 = yp - ym;
    kp = 1;
    km = 64;
    e1 = kp/km * Cm * e2_virtual;
    e2 = phi' * zeta - Cm * e2_virtual;
    eps1 = e1 + psi * e2;
    % eps1 = kp/km * phi' * zeta + psi * e2;
    u = phi' * w_bar + psi * r(t);

    % Derivative of states :
    dx = zeros(nStates, 1);

    dx(xp_idx) = Ap * xp + Bp * u;
    dx(xm_idx) = Am * xm + Bm * r(t);

    dx(w1_idx) = LAMBDA * w1 + [zeros(n-2, 1); l] * u;
    dx(w2_idx) = LAMBDA * w2 + [zeros(n-2, 1); l] * yp;

    % Adaptation law :
    dx(phi_idx) = - gamma * eps1 * zeta;
    dx(psi_idx) = - gamma * eps1 * e2;
    dx(e2_virtual_idx) = Am * e2_virtual + Bm * phi' * w_bar;

    dzeta_virtual = zeros(m, numel(w_bar));
    for i = 1:numel(w_bar)
        dzeta_virtual(:, i) = Am * zeta_virtual(:, i) + Bm * w_bar(i);
    end
    dx(zeta_virtual_idx) = dzeta_virtual(:);

end