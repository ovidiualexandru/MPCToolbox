%INV_PENDULUM_MPT Run the inverted pendulum simluation with an MPT MPC
%formulation.
clear
addpath('./pendulum');
addpath('./util');
%% System
A = [ 1.0259     0.5040   0       0;
      1.0389     1.0259   0       0;
     -0.0006          0   1    0.05;
     -0.0247    -0.0006   0       1];
B = [-0.0013; -0.0504; 0.0006; 0.025];
x0 = [0.15; 0; 0; 0]; %Initial state
N = 150; % simulation steps
Nc = 20; % control and prediction horizon
d = zeros(4,N); %disturbance vector
dist_k = 60;
d(2, dist_k) = 0.0;
%% Cost matrices and constraints
Q = [ 1      0   0     0; 
      0   0.01   0     0; 
      0      0   1     0; 
      0      0   0  0.01];
R = 0.01;
dx = [0.2, inf, 1, inf;
      -0.2, -inf, -1, -inf]; %state constraints, positive and negative
du = [5; -5]; %input constraints
%% Solver initialization
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
FVAL = zeros(1, N); %save cost value
TEVAL = zeros(1, N); %save calculation time
x = x0;
xr = x0; % 'real' x
%% MPC Solve
sys = LTISystem('A', A, 'B', B, 'Ts', 0.05); % check this?
ctrl = MPCController(sys, Nc);
ctrl.model.x.min = dx(2,:)';
ctrl.model.x.max = dx(1,:)';
ctrl.model.u.min = du(2,:)';
ctrl.model.u.max = du(1,:)';
ctrl.model.x.penalty = QuadFunction(Q);
ctrl.model.u.penalty = QuadFunction(R);
for i = 1:N
    tic;
    %% Get next command 
    u = ctrl.evaluate(x);
    teval = toc;
    %% Data logging
    X(:,i) = x; %save current state
    U(:,i) = u;
    TEVAL(i) = teval;
    %% Send to plant
    xr = A*xr + B*u + d(:,i);
    x = xr + 0.00*rand(nx,1) + 0.00*rand(nx,1).*xr;
end
%% Plotting
t = 0:N-1;
figure;
clf;
subplot(4,2,1);
plot(t,U);
rescaleYLim(gca, [du(2) du(1)]*1.1);
grid on
title('Input u');
line([0;N],[du(1);du(1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[du(2);du(2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,2);
plot(t,X(1,:));
rescaleYLim(gca, [dx(2,1) dx(1,1)]*1.1);
grid on
title('Arm position x_1');
line([0;N],[dx(1,1);dx(1,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,1);dx(2,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,4);
plot(t,X(2,:));
rescaleYLim(gca, [dx(2,2) dx(1,2)]*1.1);
grid on
title('Arm speed x_2');
line([0;N],[dx(1,2);dx(1,2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,2);dx(2,2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,6);
plot(t,X(3,:));
grid on
rescaleYLim(gca, [dx(2,3) dx(1,3)]*1.1); 
title('Trolley position x_3');
line([0;N],[dx(1,3);dx(1,3)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,3);dx(2,3)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,8);
plot(t,X(4,:));
rescaleYLim(gca, [dx(2,4) dx(1,4)]*1.1);
grid on
title('Trolley speed x_4');
line([0;N],[dx(1,4);dx(1,4)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,4);dx(2,4)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

axes('Position',[0.25 0.9 0.5 0.1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
titlestring = sprintf('Run with N=%d, Nc = %d',N,Nc);
text(0.5, 1, ['\bf ' titlestring],'HorizontalAlignment','center','VerticalAlignment', 'top')

plot_ft(FVAL, TEVAL, 'LMPC(MPT) Inverted Pendulum Performance',figure);
