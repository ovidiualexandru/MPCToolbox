%QUANSER_MPC_CONDENSED Run the Quanser 3-DOF helicopter simulation with a
%condensed MPC formulation.
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
Q = diag([5, .1, .5, .1, .1, .1],0);
R = diag([.01, .01],0);
%state constraints, positive and negative
dx = [ 30,  50,  90,  50,  inf,  inf;
      -30, -50, -90, -50, -inf, -inf];
%input constraints
du = [ 22,  22;
      -22, -22];
%% Model generation
% Set model coefficients. Leave empty for default value
mpc_param= []; % Use nominal model

sim_param.Jepsilon = []; %Default value: 0.86 kg*m^2
sim_param.Jtheta = []; %Default value: 0.044 kg*m^2
sim_param.Jphi = []; %Default value: 0.82 kg*m^2
sim_param.La = []; %Default value: 0.62 m
sim_param.Lc = []; %Default value: 0.44 m
sim_param.Ld = []; %Default value: 0.05 m
sim_param.Le = []; %Default value: 0.02 m
sim_param.Lh = []; %Default value: 0.177 m
sim_param.Mf = []; %Default value: 0.69 kg
sim_param.Mb = []; %Default value: 0.69 kg
sim_param.Mc = []; %Default value: 1.69 kg
sim_param.Km = []; %Default value: 0.5 N/V
sim_param.niu_epsilon = []; %Default value: 0.001 kg*m^2/s
sim_param.niu_theta = []; %Default value: 0.001 kg*m^2/s
sim_param.niu_phi = []; %Default value: 0.005 kg*m^2/s

%Get MPC continous model
mpc_sl = quanser_model('sl', mpc_param); 
%Get simulation continous model
sim_nl_c = quanser_model('nl', sim_param);
%Get simulation discrete model
sim_nl_d = nonlinear_c2d(sim_nl_c, h, 'euler');
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
FVAL = zeros(1, N); %save cost value
TEVAL = zeros(1, N); %save calculation time
x = x0;
xr = x0; % 'real' x
u = u0;
%% MPC solve
ue = []; %input estimated solution
for i = 1:N
    %% Update SL Model
    tic;
    if mod(i,L) == 0 || i == 1
        % Compute (A,B,g) using the MPC model
        [A,B,g] = mpc_sl(x,u); %recalculate (A,B,g)
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
    [ue, Xe,fval,EXITFLAG, OUTPUT] = lmpc_condensed(...
        Ad, Bd, Q, R, Nc, du_bar, dx_bar, xbar, xrefbar, urefbar, ue);
    if EXITFLAG < 0
        fprintf('Iteration %d\n',i)
        error('Solver error \n');
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
    xr = sim_nl_d(xr,u);
    x = xr + 0.0*rand(nx,1) + 0.0*rand(nx,1).*xr;
end
fprintf('\n');
%% Plotting
quanser_plot(X,U,dx, du,'MPC-SL(condensed form) Quanser Plot',7, XREF);
quanser_phase_plot(X, 'MPC-SL(condensed form) Quanser Phase-Plot',8, XREF);
plot_ft(FVAL, TEVAL, 'MPC-SL(condensed form) Quanser Performance',9);
