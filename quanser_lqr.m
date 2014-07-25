%QUANSER_LQR Run the Quanser 3-DOF helicopter simulation with a LQR.
clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [20; 0; 10; 0; 15; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 500; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 3; % control and prediction horizon
%% Cost matrices
Q = diag([1, .01, 1, .01, 1, .01],0);
R = diag([1, 1],0);
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
    if mod(i,Np) == 0 || i == 1
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        K = lqr(A,B,Q,R,0);
        [x_o, u_o] = affine_eq(A,B,g);
        fprintf('%d ', i);
        if mod(i,20*Np) == 0
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
    xr = quanser_disc_nl(xr,u,h);
    x = xr + 0.0*rand(nx,1) + 0.0*rand(nx,1).*xr;
end
%% Plotting
dx = [inf, inf, inf, inf, inf, inf;
      -inf, -inf, -inf, -inf, -inf, -inf]; %state constraints, positive and negative
du = [inf, inf;
      -inf, -inf]; %input constraints
quanser_plot(X,U,dx,du,'LQR Quanser Plot',4);
quanser_phase_plot(X, 'LQR Quanser Phase-Plot',5);
plot_ft(FVAL, TEVAL, 'LQR Quanser Performance', 6);
