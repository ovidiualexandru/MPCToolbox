clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [15; 0; 0; 0; 20; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 200; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 3; % control and prediction horizon
Nc = 3;
%% Cost matrices and constraints
Q = diag([1, 1, 1, 1, 1, 1],0);
R = diag([1, 1],0);
dx = [inf, inf, inf, inf, inf, inf;
      -inf, -inf, -inf, -inf, -inf, -inf]; %state constraints, positive and negative
du = [inf, inf;
      -inf, -inf]; %input constraints
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
x = x0;
u = u0;
%% MPC solve
for i = 1:N
    %% Iteration printing
    if mod(i,Np) == 0
        fprintf('%d ', i);
        if mod(i,20*Np) == 0
            fprintf('\n');
        end
    end
    %% Get next command
    [ue, Xe,FVAL,EXITFLAG, OUTPUT] = nmpc_fullspace(@quanser_disc_nl, h, Q, R, Nc, du, dx, x);
    if EXITFLAG < 0
        fprintf('Iteration: %d, EXITFLAG: %d\n',i, EXITFLAG)
        error('Solver error');
    end
    u = ue(:,1); %use only the first command in the sequence
    %% Data logging
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    %% Send to plant
    xr = quanser_disc_nl(x,u,h);
    % x = xr + 0.1.*rand(nx,1).*xr;
    x = xr;
end
%% Plotting
quanser_plot(X,U,dx, du,'Nonlinear-MPC Quanser Plot',1);
quanser_phase_plot(X, 'Nonlinear-MPC Quanser Phase-Plot',2);
%% Clean-up
rmpath('./quanser');
rmpath('./util');