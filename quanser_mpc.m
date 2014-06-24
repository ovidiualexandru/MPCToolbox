clear
%% System
x0 = [45; 0; 5; 0; 30; 0]; %Initial state
N = 30; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Nc = 5; % control and prediction horizon
%% Cost matrices and constraints
Q = diag([1, 0.1, 1, 0.1, 1, 0.1],0);
R = diag([0.2,0.2],0);
dx = [inf; inf; inf; inf; inf; inf;
      -inf; -inf; -inf; -inf; -inf; -inf]; %state constraints, positive and negative
du = [inf; inf;
     -inf; -inf]; %input constraints
%% QP solve
X = zeros(nx, N); %save all states, for plotting
u = zeros(nu, N); %save all inputs
x = x0;
C = [1 0 0 0 0 0; 0 0 0 0 1 0];
% [A,B,g] = quanser_cont_sl(x,u(:,1)); %Initial (A,B) pair
for i = 1:N
    X(:,i) = x; %save current state
    if mod(i,l) == 0
        [A,B,g] = quanser_cont_sl(x,u(:,i)); %recalculate (A,B,g)
    end
    sys = ss(A,B,C,0);
    sysd = c2d(sys,h, 'zoh');
    Ad = sysd.a;
    Bd = sysd.b;
    [ue, Xe] = qp_fullstate(Ad, Bd, Q, R, Nc, du, dx, x);
    u(:,i) = ue(:,1); %use only the first command from predictions
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u(:,i)]);
    x = Yout(end, 1:6)'; %get new state
    x = x + 0.0.*rand(nx,1).*x ; %add noise and disturbance
end
%% Plotting
t = 1:N;

figure(1);
clf;
whitebg([0 0 0]);

%Plot the input 3 times, for each state pair
for i = 1:3
    subplot(3,3,i);
    plot(t, u(1,:) ,'y--', t, u(2,:), 'c--');
    title('Inputs');
    grid on
    xlabel('[k]');
    ylabel('[volts]');
    if i == 1
        legend('Vf', 'Vb', 'Location', 'Best');
    end
end

%Plot the states
subplot(3,3,1+3);
plot(t,X(1,:), 'b-');
title('Elevation angle $\epsilon$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on
legend('NL ode45', 'NL euler', 'SL c2d', 'SL euler', 'Location', 'Best');

subplot(3,3,4+3);
plot(t,X(2,:), 'b-');
title('Elevation speed $\dot{\epsilon}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on

subplot(3,3,2+3);
plot(t,X(3,:), 'b-');
title('Pitch angle $\theta$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on

subplot(3,3,5+3);
plot(t,X(4,:), 'b-');
title('Pitch speed $\dot{\theta}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on

subplot(3,3,3+3);
plot(t,X(5,:), 'b-');
title('Travel angle $\phi$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on

subplot(3,3,6+3);
plot(t,X(6,:), 'b-');
title('Travel speed $\dot{\phi}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on