clear
%% Initialization
x0 = [45; 0; 5; 0; 30; 0]; %Initial state
N = 50; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
%% Input signal shape
u = ones(nu, N);
u(:, 1:20) = repmat([5; 5],1,20);
u(:, 21: 30) = repmat([3.7; 2.5],1,10);
u(:, 31:50) = repmat([3.5; 3.2],1,20);
%% Nonlinear model with ode45
Xtil = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xtil(:,i) = x; %save current state
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u(:,i)]);
    x = Yout(end, 1:6)'; %get new state
end
%% Nonlinear model with euler discretization
Xhat = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xhat(:,i) = x; %save current state
    f = quanser_cont_nl([],[x; u(:,i)]);
    xd = f(1:6);
    x = x + h*xd;
end
%% Succesive Liniarization discrete model
Xbard = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B] = quanser_cont_sl(x,u(:,1)); %Initial (A,B) pair
C = [1 0 0 0 0 0; 0 0 0 0 1 0];
for i = 1:N
    Xbard(:,i) = x; %save current state
    [A,B] = quanser_cont_sl(x,u(:,i)); %recalculate (A,B) for each timestep    
    sys = ss(A,B,C,0);
    sysd = c2d(sys,h);
    Ad = sysd.a;
    Bd = sysd.b;
    x = Ad*x + Bd*u(:,i);
end
% Xbard = zeros(nx, N);
%% Succesive Liniarization model with euler discretization
Xbarc = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B] = quanser_cont_sl(x,u(:,1)); %Initial (A,B) pair
for i = 1:N
    Xbarc(:,i) = x; %save current state
    [A,B] = quanser_cont_sl(x,u(:,i)); %recalculate (A,B) for each timestep
    xd = A*x + B*u(:,i);
    x = x + h*xd;
end
% Xbarc = zeros(nx, N);
%% Plotting 
% t = 0:h:((N-1)*h);
t = 1:N;

figure(1);
clf
% set(gcf, 'Units', 'normalized');
% set(gcf, 'Position', [1 0.4 1 0.5]); % set the figure at 50% height screen 2
whitebg([0 0 0]);

subplot(2,3,1);
plot(t,Xtil(1,:), 'b-');
title('Elevation angle $\epsilon$','Interpreter','latex');
hold on
plot(t,Xhat(1,:), 'r:');
plot(t,Xbard(1,:), 'c--');
plot(t,Xbarc(1,:), 'g:');
xlabel('[k]');
ylabel('[deg]');
grid on
legend('NL ode45', 'NL euler', 'SL c2d', 'SL euler', 'Location', 'Best');
hold off

subplot(2,3,4);
plot(t,Xtil(2,:), 'b-');
title('Elevation speed $\dot{\epsilon}$','Interpreter','latex');
hold on
plot(t,Xhat(2,:), 'r:');
plot(t,Xbard(2,:), 'c--');
plot(t,Xbarc(2,:), 'g:');
xlabel('[k]');
ylabel('[deg/s]');
grid on
hold off

subplot(2,3,2);
plot(t,Xtil(3,:), 'b-');
title('Pitch angle $\theta$','Interpreter','latex');
hold on
plot(t,Xhat(3,:), 'r:');
plot(t,Xbard(3,:), 'c--');
plot(t,Xbarc(3,:), 'g:');
xlabel('[k]');
ylabel('[deg]');
grid on
hold off

subplot(2,3,5);
plot(t,Xtil(4,:), 'b-');
title('Pitch speed $\dot{\theta}$','Interpreter','latex');
hold on
plot(t,Xhat(4,:), 'r:');
plot(t,Xbard(4,:), 'c--');
plot(t,Xbarc(4,:), 'g:');
xlabel('[k]');
ylabel('[deg/s]');
grid on
hold off

subplot(2,3,3);
plot(t,Xtil(5,:), 'b-');
title('Travel angle $\phi$','Interpreter','latex');
hold on
plot(t,Xhat(5,:), 'r:');
plot(t,Xbard(5,:), 'c--');
plot(t,Xbarc(5,:), 'g:');
xlabel('[k]');
ylabel('[deg]');
grid on
hold off

subplot(2,3,6);
plot(t,Xtil(6,:), 'b-');
title('Travel speed $\dot{\phi}$','Interpreter','latex');
hold on
plot(t,Xhat(6,:), 'r:');
plot(t,Xbard(6,:), 'c--');
plot(t,Xbarc(6,:), 'g:');
xlabel('[k]');
ylabel('[deg/s]');
grid on
hold off

figure(2);
clf
whitebg([0 0 0]);
plot(t, u(1,:) ,'y-', t, u(2,:), 'c--');
title('Inputs');
grid on
legend('Vf', 'Vb', 'Location', 'Best');
xlabel('[k]');
ylabel('[volts]');
