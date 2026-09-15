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
params.u0 = 0.5;   % example drug infusion rate
params.gamma = 0.9;    % example drug decay rate 
%Estrogen related parameters p4
params.lambda1= 0.20;     % Tumor formation rate as a result of DNA damage by excess estrogen
params.lambda3= 0.002;    %  Estrogen-induced immune suppression
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
y0=[5; 4; 3; 5; 0];
% Time span and step size
t0 = 0;
tf = 10;
dt = .001;
t = t0:dt:tf;
y= zeros(length(y0), length(t));
y(:, 1) = y0;
% RK4 integration
for i = 1:length(t)-1
    k1 = dt * ref_model(t(i), y(:,i), params );
    k2 = dt * ref_model(t(i)+dt/2, y(:,i)+k1/2, params);
    k3 = dt * ref_model(t(i)+dt/2, y(:,i)+k2/2, params);
    k4 = dt * ref_model(t(i)+dt, y(:,i)+k3, params);
    y(:,i+1) = y(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
t = t0:dt:tf;
v0=[5; 4; 3; 5; 0; 1; 1];
v = zeros(length(v0), length(t));
v(:, 1) = v0;
% RK4 integration
for i = 1:length(t)-1
    k1 = dt * modified_modelES(t(i), v(:,i), params);
    k2 = dt * modified_modelES(t(i)+dt/2, v(:,i)+k1/2, params);
    k3 = dt * modified_modelES(t(i)+dt/2, v(:,i)+k2/2, params);
    k4 = dt * modified_modelES(t(i)+dt, v(:,i)+k3, params);
    v(:,i+1) = v(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
figure;
plot(t,y(1, :),'c','Linewidth',6.0)
set(gca,'FontSize',16)
xlabel('\bf Time(days)','fontsize',30,'linewidth',30);ylabel('\bf  P_T Population','fontsize',30,'linewidth',30);
set(gca,'linewidth',2.5); box off;
set(gca,'FontSize',30)
hold on
plot(t, v(1, :),'r','Linewidth',6.0)
set(gca,'FontSize',16)
xlabel('\bf  Time(days) ','fontsize',30,'linewidth',30);ylabel('\bf Tumor Cell (P_T) ','fontsize',30,'linewidth',30);
set(gca,'linewidth',2.5); box off;
set(gca,'FontSize',30)
[~, hobj, ~, ~] = legend({'\bf P_T without Est. and \beta-cat. ','\bf P_T with Est. and \beta-cat.'},'Fontsize',10,'Location','northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',10);
ht = findobj(hobj,'type','text');
set(ht,'FontSize',16);
set(gca,'linewidth',2.5,'FontWeight','Bold'); box off;
set(gca,'FontSize',16)
legend boxoff       
