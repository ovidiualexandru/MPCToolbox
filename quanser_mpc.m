clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [5; 0; 0; 0; 3; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 1000; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 5; % control and prediction horizon
Nc = 5;
%% Cost matrices and constraints
Q = diag([10, .01, 10, .01, 10, .01],0);
R = diag([.01, .01],0);
dx = [90, inf, 90, inf, 90, inf;
      -90, -inf, -90, -inf, -90, -inf]; %state constraints, positive and negative
du = [5, 5;
      -5, -5]; %input constraints
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
x = x0;
u = u0;
[A,B,g] = quanser_cont_sl(x,u); %Initial (A,B,g) pair
Ad = eye(nx) + h*A;
Bd = h*B;
n = linsolve([A B], -g);
x_o = n(1:nx);
u_o = n(nx + 1:end);
du_bar = du - repmat(u_o',2,1);
dx_bar = dx - repmat(x_o',2,1);
%% MPC solve
for i = 1:N
    %% Update SL Model
    if mod(i,Np) == 0
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        Ad = eye(nx) + h*A;
        Bd = h*B;
        n = linsolve([A B], -g);
        x_o = n(1:nx);
        u_o = n(nx + 1:end);
        du_bar = du - repmat(u_o',2,1);
        dx_bar = dx - repmat(x_o',2,1);
        fprintf('%d ', i);
        if mod(i,20*Np) == 0
            fprintf('\n');
        end
    end
    %% Get next command
    xbar = x - x_o;
    [ue, Xe,FVAL,EXITFLAG, OUTPUT] = lmpc_sparse(Ad, Bd, Q, R, Nc, du_bar, dx_bar, xbar);
    if EXITFLAG < 0
        fprintf('Iteration %d\n',i)
        error('Quadprog error ');
    end
    ubar = ue(:,1); %use only the first command in the sequence
    u = ubar + u_o;
    %% Data logging
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    %% Send to plant
    xr = quanser_disc_nl(x,u,h);
    % x = xr + 0.1.*rand(nx,1).*xr;
    x = xr;
end
%% Plotting
quanser_plot(X,U,dx, du,'MPC Quanser Plot',1);
quanser_phase_plot(X, 'MPC Quanser Phase-Plot',2);
%% Clean-up
rmpath('./quanser');
rmpath('./util');