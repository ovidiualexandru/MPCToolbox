clear
%% System
x0 = [10; 0; 20; 0; 10; 0]; %Initial state
Vf = 1.8865;
Vb = 1.936;
N = 2000; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Nc = 6; % control and prediction horizon
%% Cost matrices and constraints
Q = diag([3, 10, 1, 5, 10, 2],0);
% Q = eye(6)*1;
R = diag([0.1, 0.1],0);
% R = eye(2)*1;
dx = [inf; inf; inf; inf; inf; inf;
      -inf; -inf; -inf; -inf; -inf; -inf]; %state constraints, positive and negative
du = [inf; inf;
     -inf; -inf]; %input constraints
%% LQR solve
X = zeros(nx, N); %save all states, for plotting
u = zeros(nu, N); %save all inputs
u(:,1) = [Vf; Vb];
x = x0;
C = [1 0 0 0 0 0; 0 0 0 0 1 0];
% [A,B,g] = quanser_cont_sl(x,u(:,1)); %Initial (A,B) pair
[A,B,g] = quanser_cont_sl(x,u(:,1)); %recalculate (A,B,g)
K = lqr(A,B,Q,R,0);
% u_o = u(:,1);
% x_o = A\(-g-B*u(:,1));
% x_o = linsolve(A, -g - B*u(:,1));
x_o = x;
% u_o = B\(-g - A*x);
u_o = linsolve(B, -g - A*x);
for i = 1:N
    X(:,i) = x; %save current state
    if mod(i,Nc) == 0
        [A,B,g] = quanser_cont_sl(x,u(:,i)); %recalculate (A,B,g)
        K = lqr(A,B,Q,R,0);
        %% Solution one: fix u_o
%         u_o = u(:,i);
%         x_o = A\(-g-B*u(:,i));
%         x_o = linsolve(a, -g - B*u(:,i);
        %% Solution two: fix x_o
        x_o = x;
         % u_o = B\(-g - A*x);
        u_o = linsolve(B, -g - A*x);
        fprintf('%d ', i);
        if mod(i,20*Nc) == 0
            fprintf('\n');
        end
    end
    ubar = -K*(x - x_o);
    u(:,i) = ubar + u_o;
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u(:,i)]);
    x = Yout(end, 1:6)'; %get new state
%     x = x + 0.0.*rand(nx,1).*x ; %add noise and disturbance
end
%% Plotting
t = 1:N;

figure(1);
clf;
whitebg([0 0 0]);

%Plot the input 3 times, for each state pair
for i = 1:3
    subplot(3,3,i);
    plot(t, u(1,:) ,'y--', t, u(2,:), 'c:');
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