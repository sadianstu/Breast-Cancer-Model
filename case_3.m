clc;
clear;

params.r1 = 0.5045;
params.r2 = 0.9;
params.r3 = 0.6169;
params.b1=2e-9;
params.b2=3.6e-6;
params.b3=1.83e-9;

params.k1 = 5e8;
params.k2 = 2.777e5;
params.k3 = 5.4641e8;

params.beta = 0.0937;
params.rho1 = 0.0122;

params.alpha1 = 0.014;
params.delta1 = 0.0155;
params.alpha2 = 4.34e-14;

params.sigma1 = 0.0937;
params.sigma2 = 0.0003;
params.rho2 = 4.388e-14;
params.delta2 = 0.00005;

params.alpha3 = 9e-4;
params.sigma3 = 6e-4;
params.rho3 = 0.0088;
params.delta3 = 9e-4;

params.rho = 0.05;
params.eta = 0.3;

params.u0 = 0;   % example drug infusion rate
params.gamma = 0.9;    % example drug decay rate
params.lambda1= 0.20;     % Tumor induction by estrogen
params.lambda3= 0.002;    % Estrogen suppression of immune cells
params.epsilon= 0.8;   % Estrogen input rate (pg/mL/day)
params.mu4= 0.97;    % Natural decay rate of estrogen
params.g= 15;     % Threshold for immune suppression
params.xi= 0.175;    % Drug suppression sensitivity (estimated to match P4's k = 0.6)

% Beta-catenin parameters from Paper P3
params.kS  = 30;           % ?-catenin production rate
params.dS  = 15;           % ?-catenin degradation rate
params.lambda_ST = 0.75;     % Tumor stimulation rate due to ?-catenin
params.KS  = 25;            % Half-saturation constant
params.lambda_SN = 0.2;   % New parameter for normal cell loss due to ?-catenin(estimated)
params.kh = 0.005;         % beta-catenin suppression of immune activation (estimated)

%% Control-related
params.u1 = 0; params.u2 = 0; params.u3 = 0.4;
T0=5;
H0=4;
R0=3;
N0=5;
D0=1;
E0=1;
S0=1;

% Time parameters
t0 = 0;
tf = 5;
dt = .01;
y0 = [T0; H0; R0; N0; D0; E0; S0];

% Runge-Kutta (RK4) method
t = t0:dt:tf;
y0 = [T0; H0; R0; N0; D0; E0; S0];
y = zeros(length(y0), length(t));
y(:, 1) = y0;
% RK4 integration
for i = 1:length(t)-1
    k1 = dt * model(t(i), y(:,i), params);
    k2 = dt * model(t(i)+dt/2, y(:,i)+k1/2, params);
    k3 = dt * model(t(i)+dt/2, y(:,i)+k2/2, params);
    k4 = dt * model(t(i)+dt, y(:,i)+k3, params);
    y(:,i+1) = y(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end


% Runge-Kutta (RK4) method
t = t0:dt:tf;
z0 = [T0; H0; R0; N0; D0; E0; S0];
z = zeros(length(z0), length(t));
z(:, 1) = z0;
% RK4 integration
for i = 1:length(t)-1
    k1 = dt * model_control(t(i), z(:,i), params);
    k2 = dt * model_control(t(i)+dt/2, z(:,i)+k1/2, params);
    k3 = dt * model_control(t(i)+dt/2, z(:,i)+k2/2, params);
    k4 = dt * model_control(t(i)+dt, z(:,i)+k3, params);
    z(:,i+1) = z(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

%% ==== Plot Tumor Dynamics ====
figure;
plot(t, y(4,:),'b','LineWidth',4.5); hold on;
set(gca,'FontSize',16)
xlabel('\bf Time(days)','fontsize',30,'linewidth',30);ylabel('\bf Normal Cell (P_N)','fontsize',30,'linewidth',30);
set(gca,'linewidth',2); box off;
set(gca,'FontSize',30)
hold on
plot(t, z(4,:),'r','LineWidth',4.5);
set(gca,'FontSize',16)
xlabel('\bf  Time(days) ','fontsize',30,'linewidth',30);ylabel('\bf Normal cell (P_N)','fontsize',30,'linewidth',30);
set(gca,'linewidth',2); box off;
set(gca,'FontSize',30)
[~, hobj, ~, ~] = legend({'\bf without control ','\bf with control u_3 '},'Fontsize',15,'Location','northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',8);
ht = findobj(hobj,'type','text');
set(ht,'FontSize',20);
set(gca,'linewidth',2,'FontWeight','Bold'); box off;
set(gca,'FontSize',20)
legend boxoff       
