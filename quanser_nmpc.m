%QUANSER_NMPC Run the Quanser 3-DOF helicopter simulation with a
%Nonlinear-MPC formulation
clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [5; 0; 0; 0; 0; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 10; % control and prediction horizon
Nc = 10;
%% Reference state
load('references/ref1.mat'); %load XREF and UREF into workspace
N = size(XREF,2); % Simulation size
%% Cost matrices and constraints
Q = diag([1, .1, .5, .1, .1, .1],0);
R = diag([.01, .01],0);
dx = [30, 50, 90, 50, inf, inf;
      -30, -50, -90, -50, -inf, -inf]; %state constraints, positive and negative
du = [22, 22;
      -22, -22]; %input constraints
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
FVAL = zeros(1, N); %save cost value
TEVAL = zeros(1, N); %save calculation time
x = x0;
xr = x0; % 'real' x
u = u0;
%% MPC solve
for i = 1:N
    %% Iteration printing
    tic;
    if mod(i,Np) == 0
        fprintf('%d ', i);
        if mod(i,20*Np) == 0
            fprintf('\n');
        end
    end
    %% Get next command
    idif = Nc - 1;
    if i + Nc > N
        idif = N - i;
    end
    uref = UREF(:,i:i+idif);
    xref = XREF(:,i:i+idif); % Get only as many samples as possible
    [ue, Xe,fval,EXITFLAG, OUTPUT] = nmpc_fullspace(@quanser_disc_nl_euler, h, Q, R, Nc, du, dx, x, xref, uref);
    if EXITFLAG < 0
        fprintf('Iteration: %d, EXITFLAG: %d\n',i, EXITFLAG)
        error('Solver error');
    end
    u = ue(:,1); %use only the first command in the sequence
    teval = toc;
    %% Data logging
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    FVAL(i) = fval;
    TEVAL(i) = teval;
    %% Send to plant
    xr = quanser_disc_nl(xr,u,h);
    x = xr + 0.0*rand(nx,1) + 0.0*rand(nx,1).*xr;
end
%% Plotting
quanser_plot(X,U,dx, du,'Nonlinear-MPC Quanser Plot',13, XREF);
quanser_phase_plot(X, 'Nonlinear-MPC Quanser Phase-Plot',14, XREF);
plot_ft(FVAL, TEVAL, 'Nonlinear-MPC Quanser Performance',15);
%% Trajectory save
clear XREF UREF
XREF = X;
UREF = U;
save('reference_trajectory/traj1.mat','XREF','UREF','-v7');
