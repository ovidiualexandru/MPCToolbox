%QUANSER_NMPC Run the Quanser 3-DOF helicopter simulation with a
%Nonlinear-MPC formulation
clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [0; 0; 0; 0; 0; 0]; %Initial state
u0 = [1.8; 1.8]; % [Vf Vb] initial inputs
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
L = 1; % Simulation progress update rate
Nc = 10; % Control and prediction horizon
cx = 0.0; % Internal disturbance variance
ca = 0.0; % Measurement additional white noise variance
cm = 0.0; % Measurement multiplicative white noise variance
title = 'Nonlinear-MPC';
%% Reference state
savefilename = 'simulations/sim2nl.mat';
loadfilename = 'references/ref2.mat';
trajfilename = 'references/traj2.mat';
load(loadfilename); %load XREF and UREF into workspace
% If the file is a 'path' file (not a trajectory file), set the path as a
% reference
if ~exist('XREF','var')
    XREF = XPATH;
    UREF = UPATH;
end
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
mpc_nl_c = quanser_model('nl', mpc_param); 
%Get MPC discrete model
mpc_nl_d = nonlinear_c2d(mpc_nl_c, h, 'euler');
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
Xe = []; %state estimated solution
for i = 1:N
    %% Iteration printing
    tic;
    if mod(i,L) == 0
        fprintf('%d ', i);
        if mod(i,20*L) == 0
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
    [ue, Xe,fval,EXITFLAG, OUTPUT] = nmpc_fullspace(...
        mpc_nl_d, Q, R, Nc, du, dx, x, xref, uref, Xe,ue);
    if EXITFLAG < 0
        fprintf('Iteration: %d, EXITFLAG: %d\n',i, EXITFLAG)
        error('Solver error \n');
    end
    u = ue(:,1); %use only the first command in the sequence
    teval = toc;
    %% Data logging
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    FVAL(i) = fval;
    TEVAL(i) = teval;
    %% Send to plant
    xr = sim_nl_d(xr,u) + cx*rand(nx,1);
    x = xr + ca*rand(nx,1) + cm*rand(nx,1).*xr;
end
fprintf('\n');
%% Plotting
quanser_plot(X,U,dx, du,[title ' Quanser Plot'], 13, XREF);
quanser_phase_plot(X, [title ' Quanser Phase-Plot'], 14, XREF);
plot_ft(FVAL, TEVAL, [title ' Quanser Performance'], 15);
%% Save data
simout = struct;
simout.X = X;
simout.U = U;
simout.XREF = XREF;
simout.UREF = UREF;
simout.XPATH = XPATH;
simout.UPATH = UPATH;
simout.h = h;
simout.L = L;
simout.Nc = Nc;
simout.Q = Q;
simout.R = R;
simout.dx = dx;
simout.du = du;
simout.cx = cx;
simout.ca = ca;
simout.cm = cm;
simout.mpcparam = mpc_param;
simout.simparam = sim_param;
simout.date = datestr(now);
simout.notes = title;
save(savefilename, 'simout', '-v7');
%% Trajectory save
clear XREF UREF
XREF = X;
UREF = U;
save(trajfilename, 'XREF', 'UREF', 'XPATH', 'UPATH', '-v7');
