clear
%% System initialization
x0 = [45; 2; 25; 0; 30; 0]; %Initial state
u0 = [1.8865; 1.936]; % [Vf Vb] initial inputs
N = 2000; % samples
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
Nc = 10; % control and prediction horizon
%% Cost matrices
Q = diag([10, 0.001, 10, 0.001, 10, 0.001],0);
R = diag([0.01, 0.01],0);
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
u = u0;
x = x0;
X(:,1) = x;
U(:,1) = u;
[A,B,g] = quanser_cont_sl(x,u); %Initial (A,B,g) pair
K = lqr(A,B,Q,R,0);
x_o = x;
u_o = linsolve(B, -g - A*x);
%% LQR solve
for i = 2:N
    [Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u]); %f(xk, uk)
    x = Yout(end, 1:6)'; %get new state, i.e. x = x(k)
    if mod(i,Nc) == 0
        [A,B,g] = quanser_cont_sl(x,u); %recalculate (A,B,g)
        K = lqr(A,B,Q,R,0);
        x_o = x;
        u_o = linsolve(B, -g - A*x);
        fprintf('%d ', i);
        if mod(i,20*Nc) == 0
            fprintf('\n');
        end
    end
    ubar = -K*(x - x_o);
    u = ubar + u_o; % new input
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
end
%% Plotting
tk = 1:N;
t = ((1:N)-1) * h;

figure(1);
clf;
whitebg([0 0 0]);

%Plot the input 3 times, for each state pair
for i = 1:3
    haxes = subplot(3,3,i);
    %Plot the secondary X Axis first - because of grid
    plot(t, U(1,:) ,'y--', t, U(2,:), 'c:', 'Parent', haxes);
    set(haxes, 'XAxisLocation', 'top');
    xlabel('t [s]');
    ylabel('[volts]');
    title('Inputs');
    grid on
    hold on
    %Now the primary X Axis
    haxes_pos = get(haxes, 'Position');
    haxes2 = axes('Position', haxes_pos, 'Color', 'none');
    plot(tk, U(1,:) ,'y--', tk, U(2,:), 'c:', 'Parent', haxes2);
    xlabel('samples [k]');
    grid on
    hold off
    if i == 1
        legend('Vf', 'Vb', 'Location', 'Best');
    end
    
end

%Plot the states
subplot(3,3,1+3);
plot(tk,X(1,:), 'b-');
title('Elevation angle $\epsilon$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on
legend('NL ode45', 'NL euler', 'SL c2d', 'SL euler', 'Location', 'Best');

subplot(3,3,4+3);
plot(tk,X(2,:), 'b-');
title('Elevation speed $\dot{\epsilon}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on

subplot(3,3,2+3);
plot(tk,X(3,:), 'b-');
title('Pitch angle $\theta$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on

subplot(3,3,5+3);
plot(tk,X(4,:), 'b-');
title('Pitch speed $\dot{\theta}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on

subplot(3,3,3+3);
plot(tk,X(5,:), 'b-');
title('Travel angle $\phi$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg]');
grid on

subplot(3,3,6+3);
plot(tk,X(6,:), 'b-');
title('Travel speed $\dot{\phi}$','Interpreter','latex');
xlabel('[k]');
ylabel('[deg/s]');
grid on