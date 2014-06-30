function quanser_plot(X,U,figtitle, varargin)
%% Quanser simulation plot
% Call this for a nice plot of the results from a simulation.
% Examples: quanser_mpc, quanser_lqr
%% Configuration
N = size(X,2);
t = 1:N;
Sx = {'b-', 'r-', 'b-', 'r-', 'b-', 'r-'};
Su = {'y--', 'c:'};
titles = {'Elevation angle $\epsilon$'; 'Elevation speed $\dot{\epsilon}$';
    'Pitch angle $\theta$';'Pitch speed $\dot{\theta}$';
    'Travel angle $\phi$';'Travel speed $\dot{\phi}$'};
ylabels = {'[deg]','[deg/s]','[deg]','[deg/s]','[deg]','[deg/s]'};
%% Figure initialization
rows = 3;
cols = 3;
figure(1);
clf;
set(gcf, 'Name',figtitle);
whitebg([0 0 0]);
nx = 6;
nu = 2;
if nargin > 3
    dx = varargin{1};
    du = varargin{2};
else
    dx = repmat([inf; -inf], 1, nx);
    du = repmat([inf; -inf], 1, nu);
end
%% Plot inputs
for i = 1:cols
    %% Plot state
    subplot(rows, cols, i);
    plot(t, U(1,:) ,Su{1}, t, U(2,:), Su{2});
    %% Constraint plotting
    rescaleYLim(gca, [du(2,1) du(1,1)]*1.1);
    line([0;N],[du(1,1);du(1,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
    line([0;N],[du(2,1);du(2,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
    %% Title and labels
    xlabel('samples [k]');
    ylabel('[volts]');
    title('Inputs');
    grid on
    if i == 1
        legend('Vf', 'Vb', 'Location', 'Best');
    end
end
%% Plot states
for i = 1:3
    %% Plot state
    pos = i + cols; %position in figure
    k = 2*i - 1; %state index
    subplot(rows, cols, pos );
    plot(t, X(k,:) ,Sx{k});
    %% Constraint plotting
    rescaleYLim(gca, [dx(2,k) dx(1,k)]*1.1);
    line([0;N],[dx(1,k);dx(1,k)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
    line([0;N],[dx(2,k);dx(2,k)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
    %% Title and labels
    title(titles{k},'Interpreter','latex');
    xlabel('[k]');
    ylabel(ylabels{k});
    grid on 
    if i == 1
        legend('NL ode45', 'Location', 'Best');
    end
    %% Plot secondary state - its derivative
    subplot(rows, cols, pos + cols);
    plot(t, X(k+1,:) ,Sx{k+1});
    %% Constraint plotting
    rescaleYLim(gca, [dx(2,k+1) dx(1,k+1)]*1.1);
    line([0;N],[dx(1,k+1);dx(1,k+1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
    line([0;N],[dx(2,k+1);dx(2,k+1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
    %% Title and labels
    title(titles{k+1},'Interpreter','latex');
    xlabel('[k]');
    ylabel(ylabels{k+1});
    grid on
end