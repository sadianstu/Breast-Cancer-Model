function dydt = model_control(~, y,p)
    T = y(1); % Tumor cells
    H = y(2); % Hunting immune cells
    R = y(3); % Resting immune cells
    N = y(4); % Normal cells
    D = y(5); % Drug concentration
    E= y(6);
    S= y(7);
   dTdt = p.r1*T*(1-T/p.k1)-p.alpha1*T*H-p.alpha2*T*N-p.alpha3*T*D+(p.lambda_ST*S/(p.KS+S))*T+(p.lambda1*N*E)/(1+p.xi*D);
    %dTdt = p.r1*T*(1-T/p.k1)-p.alpha1*T*H-p.alpha2*T*N-p.alpha3*T*D+(p.lambda_ST*S/(p.KS+S))*T+(p.lambda1*N*E)/(1+p.xi*D)-(p.mu4+p.xi*D+p.u2)*E;
    dHdt =(p.beta*R*H)/(1+p.kh*S)-p.sigma1*H-p.sigma2*T*H-p.sigma3*H*D;
    dRdt =p.r2*R*(1-R/p.k2)-p.rho1*R*H-p.rho2*R+(p.rho*T*R)/(T+p.eta)-p.rho3*R*D-(p.lambda3*R*E)/((p.g + E)*(1+p.xi*D));
    dNdt = p.r3*N*(1-N/p.k3)-p.delta1*T*N-p.delta2*N-p.delta3*N*D-(p.lambda1*N*E)/(1+p.xi*D)-(p.lambda_SN*S/(p.KS+S))*N;
    dDdt = p.u1-p.gamma*D;
%     dEdt = p.epsilon*(1-p.u2)-(p.mu4+p.xi*D+p.u2+p.u3)*E;
%     dSdt = p.kS(1-p.u3)-(p.dS+p.u3)*S
dEdt = p.epsilon*(1-p.u2)-(p.mu4+p.xi*D+p.u2+p.u3)*E;
    dSdt = p.kS*(1-p.u3)-(p.dS+p.u3)*S


dydt = [dTdt; dHdt; dRdt; dNdt; dDdt; dEdt; dSdt];
end