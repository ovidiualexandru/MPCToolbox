clear
%% Initialization
x0 = [20; 0; 10; 0; 5; 0]; %Initial state
N = 300; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
%% Input signal shape
u = zeros(nu, N);
u(:,11:50) = repmat([0.5; 3],1,40);
u(:, 51: 100) = repmat([3.6; 1],1,50);
u(:, 101:150) = repmat([1.1; 1.1],1,50);
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
%% Succesive Liniarization model with euler discretization
Xbar = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B] = quanser_cont_sl(x,u(:,1)); %Initial (A,B) pair
for i = 1:N
    Xbar(:,i) = x; %save current state
    [A,B] = quanser_cont_sl(x,u(:,i)); %recalculate (A,B) for each timestep
    xd = A*x + B*u(:,i);
    x = x + h*xd;
end
% Xbar = zeros(nx, N);
%% Plotting 
t = 0:h:((N-1)*h);

figure(1);
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [1 0.4 1 0.5]); % set the figure at 50% height screen 2

subplot(2,3,1);
plot(t, Xtil(1,:), 'b-');
title('Elevation angle $\epsilon$','Interpreter','latex');
hold on
plot(t,Xhat(1,:), 'r--');
plot(t,Xbar(1,:), 'g:');
grid on
legend('NL ode45', 'NL euler', 'SL', 'Location', 'Best');
hold off

subplot(2,3,4);
plot(t, Xtil(2,:), 'b-');
title('Elevation speed $\dot{\epsilon}$','Interpreter','latex');
hold on
plot(t,Xhat(2,:), 'r--');
plot(t,Xbar(2,:), 'g:');
grid on
hold off

subplot(2,3,2);
plot(t, Xtil(3,:), 'b-');
title('Pitch angle $\theta$','Interpreter','latex');
hold on
plot(t,Xhat(3,:), 'r--');
plot(t,Xbar(3,:), 'g:');
grid on
hold off

subplot(2,3,5);
plot(t, Xtil(4,:), 'b-');
title('Pitch speed $\dot{\theta}$','Interpreter','latex');
hold on
plot(t,Xhat(4,:), 'r--');
plot(t,Xbar(4,:), 'g:');
grid on
hold off

subplot(2,3,3);
plot(t, Xtil(5,:), 'b-');
title('Travel angle $\phi$','Interpreter','latex');
hold on
plot(t,Xhat(5,:), 'r--');
plot(t,Xbar(5,:), 'g:');
grid on
hold off

subplot(2,3,6);
plot(t, Xtil(6,:), 'b-');
title('Travel speed $\dot{\phi}$','Interpreter','latex');
hold on
plot(t,Xhat(6,:), 'r--');
plot(t,Xbar(6,:), 'g:');
grid on
hold off