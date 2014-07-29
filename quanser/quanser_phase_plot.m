function quanser_phase_plot(X,varargin)
%QUANSER_PHASE_PLOT produces a nice phase plot figure for the Quanser 3-DOF
%helicopter simulation.
%
%   QUANSER_PHASE_PLOT(X, FIGTITLE, FIGNUMBER) plots the simulation results
%   in X into a figure. Produces a 2-by-3 subplot with the elevation-travel
%   angles on the first line and each (state, state derivative) pair on the
%   second line, using figure FIGNUMBER and with title FIGTITLE.
%
%   Input arguments:
%   - X: 6-by-N matrix with the state snapshots
%   - figtitle: figure title. If omitted, defaults to 'Quanser Phase-Plot'
%   - fignumber: figure number. If omitted, defaults to 2.
%   - XREF: The reference state vector/matrix. If XREF is a 6-by-1 vector,
%   the reference is taken as constant. Otherwise, XREF must be a 6-by-N
%   matrix.

%% Parameter processing
N = size(X,2);
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
XREF = NaN(6,N);
if nargin > 3
    XREF = varargin{3};
    if size(XREF,2) == 1
        XREF = repmat(XREF,1,N);
    end
end
quanser_phase_plot_g(X, figtitle, fignumber, XREF);
end

function quanser_phase_plot_g(X, figtitle, fignumber, XREF)
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
plot(X(5,:), X(1,:), 'b-', XREF(5,:), XREF(1,:), 'g:');
xlabel('Travel [deg]');
ylabel('Elevation [deg]');
title('Elevation - Travel Phase Plot','Interpreter','latex');
grid on
%% Phase-plot
for i = 1:3
    k = 2*i - 1; %state index
    subplot(rows, cols, i+cols);
    plot(X(k,:), X(k+1,:), 'b-', XREF(k,:), XREF(k+1,:), 'g:');
    xlabel('Angle [deg]');
    ylabel('Angular speed [deg/s]');
    title(titles{i},'Interpreter','latex');
    grid on
end
end
