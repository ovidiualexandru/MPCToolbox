%QUANSER_LQR Run the Quanser 3-DOF helicopter simulation with a LQR.
clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [20; 0; 10; 0; 15; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 600; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
L = 3; %Simulation and linear model update rate
%% Cost matrices
Q = diag([1, .01, 1, .01, 1, .01],0);
R = diag([.01, .01],0);
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
%% LQR solve
for i = 1:N
    %% Update SL Model
    tic;
    if mod(i,L) == 0 || i == 1
        [A,B,g] = mpc_sl(x,u); %recalculate (A,B,g)
        K = lqr(A,B,Q,R,0);
        [x_o, u_o] = affine_eq(A,B,g);
        fprintf('%d ', i);
        if mod(i,20*L) == 0
            fprintf('\n');
        end
    end
    %% Get next command
    xbar = x - x_o;
    ubar = -K*xbar;
    u = ubar + u_o; % new input
    teval = toc;
    %% Data logging
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    TEVAL(i) = teval;
    %% Send to plant
    xr = sim_nl_d(xr,u);
    x = xr + 0.0*rand(nx,1) + 0.0*rand(nx,1).*xr;
end
%% Plotting
%state constraints, positive and negative
dx = [ inf,  inf,  inf,  inf,  inf,  inf;
      -inf, -inf, -inf, -inf, -inf, -inf];
%input constraints
du = [ inf,  inf;
      -inf, -inf];
quanser_plot(X,U,dx,du,'LQR Quanser Plot',4);
quanser_phase_plot(X, 'LQR Quanser Phase-Plot',5);
plot_ft(FVAL, TEVAL, 'LQR Quanser Performance', 6);
