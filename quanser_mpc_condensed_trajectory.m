clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [5; 0; 0; 0; 0; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 500; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 3; % control and prediction horizon
Nc = 3;
%% Reference state
load('trajectory.mat'); %load XREF and UREF into workspace
%% Cost matrices and constraints
Q = diag([2, .1, 1, .1, .1, .1],0);
R = diag([.01, .01],0);
dx = [30, inf, 90, inf, inf, inf;
      -30, -inf, -90, -inf, -inf, -inf]; %state constraints, positive and negative
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
    %% Update SL Model
    tic;
    if mod(i,Np) == 0 || i == 1
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        [x_o, u_o] = affine_eq(A,B,g);
        du_bar = du - repmat(u_o',2,1);
        dx_bar = dx - repmat(x_o',2,1);
        Ad = eye(nx) + h*A;
        Bd = h*B;
        fprintf('%d ', i);
        if mod(i,20*Np) == 0
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
    xref = XREF(:,i:i+idif);
    [ue, Xe,fval,EXITFLAG, OUTPUT] = lmpc_condensed(Ad, Bd, Q, R, Nc, du_bar, dx_bar, xbar, xref, urefbar);
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
quanser_plot(X,U,dx, du,'MPC-SL(condensed form) with trajectory Quanser Plot',16, XREF);
quanser_phase_plot(X, 'MPC-SL(condensed form) with trajectory Quanser Phase-Plot',17, XREF);
plot_ft(FVAL, TEVAL, 'MPC-SL(condensed form) with trajectory Quanser Performance',18);
