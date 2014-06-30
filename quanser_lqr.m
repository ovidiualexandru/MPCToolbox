clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [15; 0; 30; 0; 20; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 1000; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 5; % control and prediction horizon
%% Cost matrices
Q = diag([1, 0.01, 0.25, 0.01, 0.01, 0.01],0);
R = diag([1, 1],0);
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
x = x0;
u = u0;
[A,B,g] = quanser_cont_sl(x,u); %Initial (A,B,g) pair
K = lqr(A,B,Q,R,0);
x_o = x;
u_o = linsolve(B, -g - A*x);
%% LQR solve
for i = 1:N
    ubar = -K*(x - x_o);
    u = ubar + u_o; % new input
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u]); %f(xk, uk)
    x = Yout(end, 1:6)'; %get new state, i.e. x = x(k)
    if mod(i,Np) == 0
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        K = lqr(A,B,Q,R,0);
        x_o = x;
        u_o = linsolve(B, -g - A*x);
        fprintf('%d ', i);
        if mod(i,20*Np) == 0
            fprintf('\n');
        end
    end
end
%% Plotting
figtitle = 'LQR Quanser Plot';
quanser_plot
figtitle = 'LQR Quanser Phase-Plot';
quanser_phase_plot
%% Clean-up
rmpath('./quanser');
rmpath('./util');