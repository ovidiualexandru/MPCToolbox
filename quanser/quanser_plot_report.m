function quanser_plot_report(simout)
%QUANSER_PLOT produces a nice figure for the Quanser 3-DOF simulation.
%
%   QUANSER_PLOT(simout) plots the simulation results from X and U in 
%   figure FIGNUMBER with title FIGTITLE and reference state in XREF. This 
%   function produces a 3-by-3 subplot with the inputs on the first line
%   and each (state, state derivative) pair on each column.
%
%   Fields of simout:
%   - X: 6-by-N matrix with the state snapshots
%   - U: 2-by-N matrix of the input snapshots
%   - DX: 2-by-6 matrix with constraints on each state. Can be omitted.
%   - DU: 2-by-2 matrix with constraints on each input. Only first
%   constraint is plotted. Can be omitted.
%   - FIGTITLE: figure title. If omitted, defaults to 'Quanser Phase-Plot'
%   - FIGNUMBER: figure number. If omitted, defaults to 1.
%   - XREF: The reference state vector/matrix. If XREF is a 6-by-1 vector,
%   the reference is taken as constant. Otherwise, XREF must be a 6-by-N
%   matrix.

%% Parameter processing
dx = simout.dx;
du = simout.du;
figtitle = simout.notes;
X = simout.X;
U = simout.U;
XREF = simout.XREF;
UREF = simout.UREF;
quanser_plot_g(X, U, dx, du, figtitle, 1, XREF);
end

function quanser_plot_g(X, U, dx, du, figtitle, fignumber, XREF)
%% Configuration
N = size(X,2);
t = 1:N;
Sx = {'b-', 'r-', 'b-', 'r-', 'b-', 'r-'};
Sxref = {'r-', 'r-', 'r-', 'r-', 'r-', 'r-'};
Su = {'b--', 'r:'};
titles = {'Elevation angle $\epsilon$'; 'Elevation speed $\dot{\epsilon}$';
    'Pitch angle $\theta$';'Pitch speed $\dot{\theta}$';
    'Travel angle $\phi$';'Travel speed $\dot{\phi}$'};
ylabels = {'[deg]','[deg/s]','[deg]','[deg/s]','[deg]','[deg/s]'};
%% Figure initialization
rows = 3;
cols = 3;
figure(fignumber);
clf;
set(gcf, 'Name',figtitle);
whitebg([1 1 1]);
%% Plot inputs
figure(fignumber);
plot(t, U(1,:) ,Su{1}, t, U(2,:), Su{2});
%% Constraint plotting
rescaleYLim(gca, [du(2,1) du(1,1)]*1.1);
%Upper bound
line([0;N],[du(1,1);du(1,1)], 'LineStyle', '--', 'Color', [1 0 0]);
%Lower bound
line([0;N],[du(2,1);du(2,1)], 'LineStyle', '--', 'Color', [1 0 0]);
%% Title and labels
xlabel('samples [k]');
ylabel('[volts]');
title('Inputs');
grid on
legend('Vf', 'Vb', 'Location', 'Best');
%% Plot states
for i = 1:3
    %% Plot state
    pos = i + cols; %position in figure
    k = 2*i - 1; %state index
    figure(fignumber+i);
    plot(t, X(k,:) ,Sx{k}, t, XREF(k,:), Sxref{k});
    %% Constraint plotting
    rescaleYLim(gca, [dx(2,k) dx(1,k)]*1.1);
    %Upper bound
    line([0;N],[dx(1,k);dx(1,k)], 'LineStyle', '--', 'Color', [1 0 0]);
    %Lower bound
    line([0;N],[dx(2,k);dx(2,k)], 'LineStyle', '--', 'Color', [1 0 0]);
    %% Title and labels
    title(titles{k},'Interpreter','latex');
    xlabel('[k]');
    ylabel(ylabels{k});
    grid on 
end
end