function dvdt = modified_modelES(~, v, p)
 T = v(1); % Tumor cells
    H = v(2); % Hunting immune cells
    R = v(3); % Resting immune cells
    N = v(4); % Normal cells
    D = v(5); % Drug concentration
    E= v(6);
    S= v(7);

    
    dTdt = p.r1*T*(1-T/p.k1)-p.alpha1*T*H-p.alpha2*T*N-p.alpha3*T*D+(p.lambda_ST*S/(p.KS+S))*T+(p.lambda1*N*E)/(1+p.xi*D);
    dHdt =(p.beta*R*H)/(1+p.kh*S)-p.sigma1*H-p.sigma2*T*H-p.sigma3*H*D;
    dRdt =p.r2*R*(1-R/p.k2)-p.rho1*R*H-p.rho2*R+(p.rho*T*R)/(T+p.eta)-p.rho3*R*D-(p.lambda3*R*E)/((p.g + E)*(1+p.xi*D));
    dNdt = p.r3*N*(1-N/p.k3)-p.delta1*T*N-p.delta2*N-p.delta3*N*D-(p.lambda1*N*E)/(1+p.xi*D)-(p.lambda_SN*S/(p.KS+S))*N;
    dDdt = p.u0-p.gamma*D;
    dEdt = p.epsilon-(p.mu4+p.xi*D)*E;
    dSdt = p.kS-p.dS*S;

% Model equations
dvdt = [dTdt; dHdt; dRdt; dNdt; dDdt; dEdt; dSdt];
end