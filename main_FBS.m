clc;
clear;
close all;

%% ============================================================
%  PARAMETERS
% =============================================================

params.r1 = 0.5045;
params.r2 = 0.9;
params.r3 = 0.6169;

params.b1 = 2e-9;
params.b2 = 3.6e-6;
params.b3 = 1.83e-9;

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

params.gamma = 0.9;
params.lambda3 = 0.002;
params.epsilon = 0.8;
params.mu4 = 0.97;
params.g = 15;
params.xi = 0.175;

params.kS = 30;
params.dS = 15;
params.lambda_ST = 0.75;
params.KS = 25;
params.lambda_SN = 0.2;
params.kh = 0.005;

% Tumor induction by estrogen
params.lambda1 = 0.20;

%% ============================================================
% INITIAL CONDITIONS
% =============================================================

T0 = 5;
H0 = 4;
R0 = 3;
N0 = 5;
D0 = 1;
E0 = 1;
S0 = 1;

Y0 = [T0; H0; R0; N0; D0; E0; S0];

%% ============================================================
% OPTIMAL CONTROL PARAMETERS
% =============================================================

params.c1 = 1;
params.c2 = 1;
params.c3 = 1;
params.c4 = 1;

% Control bounds
params.umin = 0;
params.umax = 1;

% FBS parameters
params.tol = 1e-6;
params.maxIter = 500;

% Relaxation parameter
params.relax = 0.5;

%% ============================================================
% TIME PARAMETERS
% =============================================================

t0 = 0;
tf = 10;

% Baseline step size
dt = 0.01;

%% ============================================================
% CASE 1: OPTIMIZE u1 ONLY
% =============================================================

fprintf('\n====================================================\n');
fprintf('CASE 1: OPTIMAL CONTROL u_1 ONLY\n');
fprintf('====================================================\n');

case1 = FBS_solver(t0,tf,dt,Y0,params,1);

%% ============================================================
% CASE 2: OPTIMIZE u2 ONLY
% =============================================================

fprintf('\n====================================================\n');
fprintf('CASE 2: OPTIMAL CONTROL u_2 ONLY\n');
fprintf('====================================================\n');

case2 = FBS_solver(t0,tf,dt,Y0,params,2);

%% ============================================================
% CASE 3: OPTIMIZE u3 ONLY
% =============================================================

fprintf('\n====================================================\n');
fprintf('CASE 3: OPTIMAL CONTROL u_3 ONLY\n');
fprintf('====================================================\n');

case3 = FBS_solver(t0,tf,dt,Y0,params,3);

%% ============================================================
% CASE 4: OPTIMIZE ALL CONTROLS
% =============================================================

fprintf('\n====================================================\n');
fprintf('CASE 4: OPTIMAL CONTROLS u_1,u_2,u_3\n');
fprintf('====================================================\n');

case4 = FBS_solver(t0,tf,dt,Y0,params,4);

%% ============================================================
% STEP-SIZE SENSITIVITY
% =============================================================

fprintf('\n====================================================\n');
fprintf('STEP-SIZE SENSITIVITY ANALYSIS\n');
fprintf('====================================================\n');

dt_values = [0.02 0.01 0.005];

sensitivity = zeros(length(dt_values),6);

for q = 1:length(dt_values)

    dtq = dt_values(q);

    fprintf('\nRunning dt = %.4f\n',dtq);

    result = FBS_solver(t0,tf,dtq,Y0,params,4);

    sensitivity(q,1) = dtq;
    sensitivity(q,2) = length(result.t)-1;
    sensitivity(q,3) = result.J;
    sensitivity(q,4) = result.Y(1,end);
    sensitivity(q,5) = result.Y(4,end);
    sensitivity(q,6) = result.iterations;

end

%% ============================================================
% STEP-SIZE TABLE
% =============================================================

SensitivityTable = array2table(sensitivity,...
    'VariableNames',...
    {'dt','Nsteps','Jstar','PT_final','PN_final','FBS_iterations'});

disp(' ');
disp('STEP-SIZE SENSITIVITY TABLE');
disp(SensitivityTable);

writetable(SensitivityTable,'StepSizeSensitivity.csv');

%% ============================================================
% PLOT CONVERGENCE FOR CASE 4
% =============================================================

figure;

semilogy(case4.convergence,'LineWidth',2);
xlabel('FBS Iteration',...
    'FontSize',15,...
    'FontWeight','bold',...
    'FontName','Times New Roman');
ax=gca;

ax.FontSize=12;
ax.FontName='Times New Roman';
ax.FontWeight='bold';
ax.LineWidth=2.8;
box off;
grid off;
ylabel('Relative Error',...
    'FontSize',15,...
    'FontWeight','bold',...
    'FontName','Times New Roman');

title('Convergence of the Forward-Backward Sweep Method',...
  'FontSize',8,...
    'FontWeight','bold',...
    'FontName','Times New Roman');
ax=gca;

ax.FontSize=12;
ax.FontName='Times New Roman';
ax.FontWeight='bold';
ax.LineWidth=2.8;
box off;
grid off;
ax.YMinorTick = 'off';
ax.YMinorGrid = 'off';
%% ============================================================
% PLOT OPTIMAL CONTROLS
% =============================================================


%% ============================================================
% DISPLAY FINAL RESULTS
% =============================================================

fprintf('\n====================================================\n');
fprintf('FINAL RESULTS - CASE 4\n');
fprintf('====================================================\n');

fprintf('FBS iterations = %d\n',case4.iterations);
fprintf('Final relative error = %.6e\n',case4.finalError);
fprintf('Optimal objective J* = %.8f\n',case4.J);
fprintf('Final tumor population P_T(tf) = %.8f\n',case4.Y(1,end));
fprintf('Final normal cell population P_N(tf) = %.8f\n',case4.Y(4,end));
