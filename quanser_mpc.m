clear
%% System initialization
x0 = [15; 0; 30; 0; 20; 0]; %Initial state
u0 = [2; 2]; % [Vf Vb] initial inputs
N = 1000; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 5; % control and prediction horizon
%% Cost matrices and constraints
Q = diag([1, 0.01, 1, 0.01, 1, 0.01],0);
R = diag([0.3, 0.3],0);
dx = [inf; inf; inf; inf; inf; inf;
      -inf; -inf; -inf; -inf; -inf; -inf]; %state constraints, positive and negative
du = [inf; inf;
      -inf; -inf]; %input constraints
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
x = x0;
u = u0;
[A,B,g] = quanser_cont_sl(x,u); %Initial (A,B,g) pair
C = [1 0 0 0 0 0; 0 0 0 0 1 0];
sys = ss(A,B,C,0);
sysd = c2d(sys,h, 'zoh');
sysd.a = eye(nx) + h*A;
sysd.b = h*B;
Ad = sysd.a;
Bd = sysd.b;
x_o = x;
u_o = linsolve(B, -g - A*x);
%% MPC solve
for i = 1:N
    [ue, Xe] = qp_fullstate(Ad, Bd, Q, R, Np, du, dx, x);
    ubar = ue(:,1); %use only the first command in the sequence
    u = ubar + u_o;
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u]); %f(xk, uk)
    x = Yout(end, 1:6)'; %get new state, i.e. x = x(k)
    if mod(i,Np) == 0
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        sys = ss(A,B,C,0);
        sysd = c2d(sys,h, 'zoh');
        sysd.a = eye(nx) + h*A;
        sysd.b = h*B;
        Ad = sysd.a;
        Bd = sysd.b;
        x_o = x;
        u_o = linsolve(B, -g - A*x);
%         fprintf('%d ', i);
%         if mod(i,20*Np) == 0
%             fprintf('\n');
%         end
    end
end
%% Plotting
tk = 1:N;

figure(1);
clf;
whitebg([0 0 0]);

Sx = {'b-', 'b-', 'b-', 'r-', 'r-', 'r-'};
Su = {'y--', 'c:'};
[xhandles, uhandles] =  plotsim(X([1 3 5 2 4 6] ,:),U,tk,2,Sx,Su);
% Now the states are [epsilon, theta phi, epsilon_dot, theta_dot, phi_dot]
axes(uhandles(1));
legend('Vf', 'Vb', 'Location', 'Best');

for i = 1:length(uhandles)
    axes(uhandles(i));
    xlabel('samples [k]');
    ylabel('[volts]');
    title('Inputs');
end

axes(xhandles(1));
title('Elevation angle $\epsilon$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
legend('NL ode45', 'NL euler', 'SL c2d', 'SL euler', 'Location', 'Best');

axes(xhandles(2));
title('Pitch angle $\theta$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on

axes(xhandles(3));
title('Travel angle $\phi$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on

axes(xhandles(4));
title('Elevation speed $\dot{\epsilon}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on

axes(xhandles(5));
title('Pitch speed $\dot{\theta}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on

axes(xhandles(6));
title('Travel speed $\dot{\phi}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on