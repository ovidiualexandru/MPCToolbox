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
    [Tout, Yout] = ode45(@quanser_nonlin_cont, [0 h], [x; u(:,i)]);
    x = Yout(end, 1:6)'; %get new state
end
%% Nonlinear model with euler discretization
Xhat = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xhat(:,i) = x; %save current state
    [F, G] = quanser_nonlin_disc(x);
    x = x + h*(F + G*u(:,i));
end
%% Succesive Liniarization model with euler discretization
Xbar = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B] = quanser_sl_cont(x,u(:,1)); %Initial (A,B) pair
for i = 1:N
    Xbar(:,i) = x; %save current state
    [A,B] = quanser_sl_cont(x,u(:,i)); %uncomment this to recalculate
%     (A,B) for each timestep
    x = x + h*(A*x + B*u(:,i));
end
Xbar = zeros(nx, N);
%% Plotting 
figure(1);
t = 0:h:((N-1)*h);
subplot(6,1,1);
plot(t, Xtil(1,:), 'b-');
title('Elevation angle $\epsilon$','Interpreter','latex');
hold on
plot(t,Xhat(1,:), 'r--');
plot(t,Xbar(1,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,2);
plot(t, Xtil(2,:), 'b-');
title('Elevation speed $\dot{\epsilon}$','Interpreter','latex');
hold on
plot(t,Xhat(2,:), 'r--');
plot(t,Xbar(2,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,3);
plot(t, Xtil(3,:), 'b-');
title('Pitch angle $\theta$','Interpreter','latex');
hold on
plot(t,Xhat(3,:), 'r--');
plot(t,Xbar(3,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,4);
plot(t, Xtil(4,:), 'b-');
title('Pitch speed $\dot{\theta}$','Interpreter','latex');
hold on
plot(t,Xhat(4,:), 'r--');
plot(t,Xbar(4,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,5);
plot(t, Xtil(5,:), 'b-');
title('Travel angle $\phi$','Interpreter','latex');
hold on
plot(t,Xhat(5,:), 'r--');
plot(t,Xbar(5,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,6);
plot(t, Xtil(6,:), 'b-');
title('Travel speed $\dot{\phi}$','Interpreter','latex');
hold on
plot(t,Xhat(6,:), 'r--');
plot(t,Xbar(6,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off