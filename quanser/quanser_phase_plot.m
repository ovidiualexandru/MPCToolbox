function quanser_phase_plot(X,varargin)
%QUANSER_PHASE_PLOT produces a nice phase plot figure for the Quanser 3-DOF
%helicopter simulation.
%   QUANSER_PHASE_PLOT(X, figtitle, fignumber) plot the simulation results
%   in X into a figure. Produces a 2-by-3 subplot with the elevation-travel
%   angles on the first line and each (state, state derivative) pair on the
%   second line.
%
%   Arguments:
%   - X: 6-by-N matrix with the state snapshots
%   - figtitle: figure title. If omitted, defaults to 'Quanser Phase-Plot'
%   - fignumber: figure number. If omitted, defaults to 2.

%% Parameter processing
if nargin > 1
    figtitle = varargin{1};
else
    figtitle = 'Quanser Phase Plot';
end
if nargin > 2
    fignumber = varargin{2};
else
    fignumber = 2;
end
quanser_phase_plot_g(X, figtitle, fignumber);
end

function quanser_phase_plot_g(X, figtitle, fignumber)
%% Configuration
titles = {'Elevation $\epsilon$'; 'Pitch $\theta$'; 'Travel $\phi$'};
%% Figure initialization
rows = 2;
cols = 3;
figure(fignumber);
clf;
whitebg([1 1 1]);
set(gcf, 'Name',figtitle);
%% Phase-plot
subplot(rows,cols, 1:3);
plot(X(5,:), X(1,:), 'b-');
xlabel('Travel [deg]');
ylabel('Elevation [deg]');
title('Elevation - Travel Phase Plot','Interpreter','latex');
grid on
%% Phase-plot
for i = 1:3
    k = 2*i - 1; %state index
    subplot(rows, cols, i+cols);
    plot(X(k,:), X(k+1,:), 'b-');
    xlabel('Angle [deg]');
    ylabel('Angular speed [deg/s]');
    title(titles{i},'Interpreter','latex');
    grid on
end
end
