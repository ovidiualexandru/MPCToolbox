%QUANSER_MPC_SPARSE Run the Quanser 3-DOF helicopter simulation with a
%sparse MPC formulation.
clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [5; 0; 0; 0; 0; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
L = 3; %Simulation and linear model update rate
Nc = 3; %Control and prediction horizon
%% Reference state
load('references/ref1.mat'); %load XREF and UREF into workspace
N = size(XREF,2); % Simulation size
%% Cost matrices and constraints
Q = diag([1, .1, .5, .1, .1, .1],0);
R = diag([.01, .01],0);
%state constraints, positive and negative
dx = [ 30,  50,  90,  50,  inf,  inf;
      -30, -50, -90, -50, -inf, -inf];
%input constraints
du = [ 22,  22;
      -22, -22];
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
    %% Update SL Model
    tic;
    if mod(i,L) == 0 || i == 1
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        [x_o, u_o] = affine_eq(A,B,g);
        du_bar = du - repmat(u_o',2,1);
        dx_bar = dx - repmat(x_o',2,1);
        Ad = eye(nx) + h*A;
        Bd = h*B;
        fprintf('%d ', i);
        if mod(i,20*L) == 0
            fprintf('\n');
        end
    end
    %% Get next command
    xbar = x - x_o;
    idif = Nc - 1;
    if i + Nc > N
        idif = N - i;
    end
    urefbar = UREF(:,i:i+idif) - repmat(u_o,[1 idif+1]);
    xrefbar = XREF(:,i:i+idif) - repmat(x_o,[1 idif+1]);
    [ue, Xe,fval,EXITFLAG, OUTPUT] = lmpc_sparse(...
        Ad, Bd, Q, R, Nc, du_bar, dx_bar, xbar, xrefbar, urefbar);
    if EXITFLAG < 0
        fprintf('Iteration %d\n',i)
        error('Quadprog error ');
    end
    ubar = ue(:,1); %use only the first command in the sequence
    u = ubar + u_o;
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
quanser_plot(X,U,dx, du,'MPC-SL(sparse) Quanser Plot',1, XREF);
quanser_phase_plot(X, 'MPC-SL(sparse) Quanser Phase-Plot',2, XREF);
plot_ft(FVAL, TEVAL, 'MPC-SL(sparse) Quanser Performance',3);
