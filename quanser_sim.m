clear
%% Initialization
x0 = [20; 0; 10; 0; 5; 0]; %Initial state
Vf = 1.2; %V
Vb = 3.7; %V
t = 1; %s - simulation time
%% Define inputs
% u(:,11:30) = repmat([1; 3],1,20);
% u(:, 31: 50) = repmat([3; 1],1,20);
% u(:, 51:100) = repmat([1; 1],1,50);
% u(:, 101:150) = repmat([2;2],1,50);

%% Continous nonlinear model
[Tout, Yout] = ode45(@quanser_nonlin_cont, [0 t], [x0; Vf; Vb]);
Xc = Yout(:,1:6)';
%% Discrete nonlinear model - Euler discretization
hn = 0.1; % s
Nn = round(t/hn); %simulation steps
% N = size(Xc, 2);
% h = t/N;
nx = length(x0);
nu = 2;
Xn = zeros(nx, Nn); %save all states, for plotting
u = zeros(nu, Nn)+repmat([Vf; Vb], 1,Nn);
x = x0;
for i = 1:Nn
    Xn(:,i) = x; %save current state
    x = quanser_nonlin_disc(hn,x,u(:,i)); %get next state
end
%% Discrete linear model - Successive linearizations
hl = 0.01; % s
Nl = round(t/hl); %simulation steps
% N = size(Xc, 2);
% h = t/N;
nx = length(x0);
nu = 2;
Xl = zeros(nx, Nl); %save all states, for plotting
u = zeros(nu, Nl)+repmat([Vf; Vb], 1,Nl);
x = x0;
for i = 1:Nl
    Xl(:,i) = x; %save current state
    x = quanser_sl_cont(hl,x,u(:,i)); %get next state
end
%% Plotting
figure(1);
plot(linspace(0,t,Nn),Xn)
legend('Elevation angle', 'Elevation speed', 'Pitch angle', 'Pitch speed',...
    'Travel angle', 'Travel speed');
grid on
figure(2);
plot(Tout, Xc);
legend('Elevation angle', 'Elevation speed', 'Pitch angle', 'Pitch speed',...
    'Travel angle', 'Travel speed');
grid on
figure(3);
plot(linspace(0,t,Nl),Xl)
legend('Elevation angle', 'Elevation speed', 'Pitch angle', 'Pitch speed',...
    'Travel angle', 'Travel speed');
grid on