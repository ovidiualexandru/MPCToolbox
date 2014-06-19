clear
%% Initialization
x0 = [20; 0; 10; 0; 5; 0]; %Initial state
N = 150; % samples
h = 0.05; % s - sampling time
nu = 2;
nx = length(x0);
%% Define inputs
u = zeros(nu, N);
u(:,11:30) = repmat([1; 3],1,20);
u(:, 31: 50) = repmat([3; 1],1,20);
u(:, 51:100) = repmat([1; 1],1,50);
u(:, 101:150) = repmat([2;2],1,50);

%% Continous nonlinear model
Xc = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xc(:,i) = x; %save current state
    [Tout, Yout] = ode45(@quanser_nonlin_cont, [0 h], [x; u(:,i)]);
    x = Yout(end, 1:6)'; %get next state
end
%% Discrete nonlinear model - Euler discretization
Xn = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xn(:,i) = x; %save current state
    [F, G] = quanser_nonlin_disc(x); %get next state
    x = x + h*(F + G*u(:,i));
end
%% Discrete linear model - Successive linearizations
Xl = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B] = quanser_sl_cont(x,u(:,1)); %Initial (A,B) pair
for i = 1:N
    Xl(:,i) = x; %save current state
    [A,B] = quanser_sl_cont(x,u(:,i)); %uncomment this to recalculate
%     (A,B) for each timestep
    x = x + h*(A*x + B*u(:,i));
end
%% Plotting 
figure(1);
t = 0:h:((N-1)*h);
subplot(6,1,1);
plot(t, Xc(1,:), 'b-');
hold on
plot(t,Xn(1,:), 'r--');
plot(t,Xl(1,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,2);
plot(t, Xc(2,:), 'b-');
hold on
plot(t,Xn(2,:), 'r--');
plot(t,Xl(2,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,3);
plot(t, Xc(3,:), 'b-');
hold on
plot(t,Xn(3,:), 'r--');
plot(t,Xl(3,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,4);
plot(t, Xc(4,:), 'b-');
hold on
plot(t,Xn(4,:), 'r--');
plot(t,Xl(4,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,5);
plot(t, Xc(5,:), 'b-');
hold on
plot(t,Xn(5,:), 'r--');
plot(t,Xl(5,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off

subplot(6,1,6);
plot(t, Xc(6,:), 'b-');
hold on
plot(t,Xn(6,:), 'r--');
plot(t,Xl(6,:), 'g:');
legend('NL ode45', 'NL euler', 'SL', 'Location', 'NorthEastOutside');
hold off