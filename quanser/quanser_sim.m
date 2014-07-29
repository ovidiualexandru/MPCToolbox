clear
%QUANSER_SIM is a script that simulates the Quanser 3-DOF helicopter 
%and plots the results in figure 1.
%
% This file needs the following files to be in the same directory:
%   - quanser_params.m: load model coefficients
%   - quanser_cont_nl.m: nonlinear model state derivative estimation
%   - quanser_con_sl.m: successive linearizations model
%
% Note: this file is obsolete and will not work.

%% Initialization
x0 = [45; 0; 5; 0; 30; 0]; %Initial state
N = 50; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Np = 5; %SL horizon i.e. how many steps until a new affine term is computed
%% Input signal shape
U = ones(nu, N);
% u(:, 1:30) = repmat([3; 3],1,30);
% u(:, 31: 60) = repmat([5; 0],1,30);
% u(:, 61:90) = repmat([0; 5],1,30);
U(1,1:5) = [2 3 2.3 0 3];
U(2,1:5) = [2.2 2.1 1.3 4 3];
U(:,6:10) = ones(2,5)*3;
%% Nonlinear model with ode45
Xtil = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xtil(:,i) = x; %save current state
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; U(:,i)]);
    x = Yout(end, 1:6)'; %get new state
end
%% Nonlinear model with euler discretization
Xhat = zeros(nx, N); %save all states, for plotting
x = x0;
for i = 1:N
    Xhat(:,i) = x; %save current state
    f = quanser_cont_nl([],[x; U(:,i)]);
    xd = f(1:6);
    x = x + h*xd;
end
%% Succesive Liniarization model, discretized with c2d 
Xbard = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B,g] = quanser_cont_sl(x,U(:,1)); %Initial (A,B) pair
C = [1 0 0 0 0 0; 0 0 0 0 1 0];
for i = 1:N
    Xbard(:,i) = x; %save current state
    if mod(i,Np) == 0
        [A,B,g] = quanser_cont_sl(x,U(:,i)); %recalculate (A,B,g)
    end
    sys = ss(A,B,C,0);
    sysd = c2d(sys,h, 'zoh');
    sysd.a = eye(nx) + h*A;
    sysd.b = h*B;
    Ad = sysd.a;
    Bd = sysd.b;
    x = Ad*x + Bd*U(:,i) + h*g;
end
% Xbard = zeros(nx, N);
%% Succesive Liniarization model with euler discretization
Xbarc = zeros(nx, N); %save all states, for plotting
x = x0;
[A,B,g] = quanser_cont_sl(x,U(:,1)); %Initial (A,B) pair
for i = 1:N
    Xbarc(:,i) = x; %save current state
    if mod(i,Np) == 0
        [A,B,g] = quanser_cont_sl(x,U(:,i)); %recalculate (A,B,g)
    end
    xd = A*x + B*U(:,i) + g;
    x = x + h*xd;
end
% Xbarc = zeros(nx, N);
%% Plotting
% Configuration
t = 1:N;
Sx = {'b-', 'r-', 'b-', 'r-', 'b-', 'r-'};
Su = {'y--', 'c:'};
titles = {'Elevation angle $\epsilon$'; 'Elevation speed $\dot{\epsilon}$';
    'Pitch angle $\theta$';'Pitch speed $\dot{\theta}$';
    'Travel angle $\phi$';'Travel speed $\dot{\phi}$'};
ylabels = {'[deg]','[deg/s]','[deg]','[deg/s]','[deg]','[deg/s]'};
% Figure initialization
rows = 3;
cols = 3;
figure(1);
clf;
whitebg([0 0 0]);
%Plot inputs
for i = 1:cols
    subplot(rows, cols, i);
    plot(t, U(1,:) ,Su{1}, t, U(2,:), Su{2});
    xlabel('samples [k]');
    ylabel('[volts]');
    title('Inputs');
    grid on
    if i == 1
        legend('Vf', 'Vb', 'Location', 'Best');
    end
end

% Plot states
for i = 1:3
    %Plot state
    pos = i + cols; %position in figure
    k = 2*i - 1; %state index
    subplot(rows, cols, pos );
    hold on
    plot(t,Xtil(k,:), 'b-');
    plot(t,Xhat(k,:), 'r:');
    plot(t,Xbard(k,:), 'c--');
    plot(t,Xbarc(k,:), 'g:');
    xlabel('[k]');
    ylabel('[deg]');
    grid on
    hold off
    title(titles{k},'Interpreter','latex');
    xlabel('[k]');
    ylabel(ylabels{k});
    grid on 
    if i == 1
        legend('NL ode45', 'NL euler', 'SL c2d', 'SL euler', ...
            'Location', 'Best');
    end
    %Plot secondary state - its derivative
    subplot(rows, cols, pos + cols);
    hold on
    plot(t,Xtil(k+1,:), 'b-');
    plot(t,Xhat(k+1,:), 'r:');
    plot(t,Xbard(k+1,:), 'c--');
    plot(t,Xbarc(k+1,:), 'g:');
    xlabel('[k]');
    ylabel('[deg/s]');
    grid on
    hold off
    title(titles{k+1},'Interpreter','latex');
    xlabel('[k]');
    ylabel(ylabels{k+1});
    grid on
end
