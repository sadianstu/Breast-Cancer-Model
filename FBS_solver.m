function result = FBS_solver(t0,tf,dt,Y0,params,caseNo)

%% ============================================================
% TIME GRID
% =============================================================

t = t0:dt:tf;
n = length(t);

%% ============================================================
% INITIAL CONTROL GUESS
% =============================================================

u1 = zeros(1,n);
u2 = zeros(1,n);
u3 = zeros(1,n);

switch caseNo

    case 1
        u1(:) = 0.5;

    case 2
        u2(:) = 0.5;

    case 3
        u3(:) = 0.5;

    case 4
        u1(:) = 0.5;
        u2(:) = 0.5;
        u3(:) = 0.5;

end

%% ============================================================
% INITIAL ADJOINT
% =============================================================

A = zeros(7,n);

convergence = zeros(params.maxIter,1);

fprintf('\n');
fprintf(' Iteration       Relative Error             J\n');
fprintf('----------------------------------------------------\n');

%% ============================================================
% FBS ITERATION
% =============================================================

for iter = 1:params.maxIter

    u1_old = u1;
    u2_old = u2;
    u3_old = u3;

    %% --------------------------------------------------------
    % FORWARD SWEEP
    % ---------------------------------------------------------

    Y = zeros(7,n);
    Y(:,1) = Y0;

    for j = 1:n-1

        Yj = Y(:,j);

        Uj = [u1(j);u2(j);u3(j)];

        Uhalf = [ ...
            (u1(j)+u1(j+1))/2;
            (u2(j)+u2(j+1))/2;
            (u3(j)+u3(j+1))/2 ];

        Uend = [u1(j+1);u2(j+1);u3(j+1)];

        k1 = dt * controlled_model(t(j),...
            Yj,Uj,params);

        k2 = dt * controlled_model(...
            t(j)+dt/2,...
            Yj+k1/2,...
            Uhalf,params);

        k3 = dt * controlled_model(...
            t(j)+dt/2,...
            Yj+k2/2,...
            Uhalf,params);

        k4 = dt * controlled_model(...
            t(j)+dt,...
            Yj+k3,...
            Uend,params);

        Y(:,j+1) = Yj + ...
            (k1+2*k2+2*k3+k4)/6;

    end

    %% --------------------------------------------------------
    % BACKWARD SWEEP
    % ---------------------------------------------------------

    A = zeros(7,n);

    % Terminal condition
    A(:,n) = zeros(7,1);

    for j = n:-1:2

        Aj = A(:,j);

        Yj = Y(:,j);

        Uj = [u1(j);u2(j);u3(j)];

        Yhalf = (Y(:,j)+Y(:,j-1))/2;

        Uhalf = [ ...
            (u1(j)+u1(j-1))/2;
            (u2(j)+u2(j-1))/2;
            (u3(j)+u3(j-1))/2 ];

        Uprev = [u1(j-1);u2(j-1);u3(j-1)];

        k1 = -dt * adjoint_model(...
            t(j),Yj,Aj,Uj,params);

        k2 = -dt * adjoint_model(...
            t(j)-dt/2,...
            Yhalf,...
            Aj+k1/2,...
            Uhalf,params);

        k3 = -dt * adjoint_model(...
            t(j)-dt/2,...
            Yhalf,...
            Aj+k2/2,...
            Uhalf,params);

        k4 = -dt * adjoint_model(...
            t(j)-dt,...
            Y(:,j-1),...
            Aj+k3,...
            Uprev,params);

        A(:,j-1) = Aj + ...
            (k1+2*k2+2*k3+k4)/6;

    end

    %% --------------------------------------------------------
    % CONTROL UPDATE
    % ---------------------------------------------------------

    u1_new = u1;
    u2_new = u2;
    u3_new = u3;

    %% u1

    if caseNo == 1 || caseNo == 4

        u1_new = -A(5,:)/params.c2;

        u1_new = min(params.umax,...
            max(params.umin,u1_new));

    end

    %% u2

    if caseNo == 2 || caseNo == 4

        u2_new = ...
            A(6,:).*(params.epsilon + Y(6,:)) ...
            /params.c3;

        u2_new = min(params.umax,...
            max(params.umin,u2_new));

    end

    %% u3

    if caseNo == 3 || caseNo == 4

        u3_new = ...
            (A(6,:).*Y(6,:) + ...
            A(7,:).*(params.kS + Y(7,:))) ...
            /params.c4;

        u3_new = min(params.umax,...
            max(params.umin,u3_new));

    end

    %% --------------------------------------------------------
    % RELAXATION
    % ---------------------------------------------------------

    theta = params.relax;

    u1 = theta*u1_new + (1-theta)*u1_old;

    u2 = theta*u2_new + (1-theta)*u2_old;

    u3 = theta*u3_new + (1-theta)*u3_old;

    %% --------------------------------------------------------
    % RELATIVE ERROR
    % ---------------------------------------------------------

    err1 = norm(u1-u1_old,inf) / ...
        max(norm(u1_old,inf),1e-12);

    err2 = norm(u2-u2_old,inf) / ...
        max(norm(u2_old,inf),1e-12);

    err3 = norm(u3-u3_old,inf) / ...
        max(norm(u3_old,inf),1e-12);

    err = max([err1 err2 err3]);

    convergence(iter) = err;

    %% --------------------------------------------------------
    % OBJECTIVE FUNCTION
    % ---------------------------------------------------------

    integrand = ...
        params.c1*Y(1,:) + ...
        0.5*params.c2*u1.^2 + ...
        0.5*params.c3*u2.^2 + ...
        0.5*params.c4*u3.^2;

    J = trapz(t,integrand);

    fprintf('%8d       %14.6e       %14.8f\n',...
        iter,err,J);

    %% --------------------------------------------------------
    % CONVERGENCE CHECK
    % ---------------------------------------------------------

    if err < params.tol

        fprintf('----------------------------------------------------\n');
        fprintf('Converged after %d iterations.\n',iter);
        fprintf('Final relative error = %.6e\n',err);

        break;

    end

end

%% ============================================================
% OUTPUT
% =============================================================

if iter == params.maxIter && err >= params.tol

    fprintf('----------------------------------------------------\n');
    fprintf('Maximum number of iterations reached.\n');
    fprintf('Final relative error = %.6e\n',err);

end

convergence = convergence(1:iter);

result.t = t;
result.Y = Y;
result.A = A;

result.U = [u1;u2;u3];

result.J = J;

result.iterations = iter;
result.finalError = err;

result.convergence = convergence;

end