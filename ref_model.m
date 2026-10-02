function dydt = ref_model(~, y,p)
    T = y(1); % Tumor cells
    H = y(2); % Hunting immune cells
    R = y(3); % Resting immune cells
    N = y(4); % Normal cells
    D = y(5); % Drug concentration

    dTdt = p.r1 * T * (1 - T/p.k1) - p.alpha1*T*H - p.alpha2*T*N - p.alpha3*T*D;
    dHdt = p.beta*R*H - p.sigma1*H - p.sigma2*T*H - p.sigma3*H*D;
    dRdt = p.r2 * R * (1 - R/p.k2) - p.rho1*R*H - p.rho2*R + (p.rho*T*R)/(T + p.eta) - p.rho3*R*D;
    dNdt = p.r3 * N * (1 - N/p.k3) - p.delta1*T*N - p.delta2*N - p.delta1*N*D;
    dDdt = p.u0 - p.gamma*D;


    dydt = [dTdt; dHdt; dRdt; dNdt; dDdt];
end