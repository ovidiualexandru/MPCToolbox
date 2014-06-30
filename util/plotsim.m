function [xhandles, uhandles]= plotsim(X, U, t, rows, Sx, Su)
%Plot a nice figure for the inputs and states
% X - matrix with states, each column a state snapshot
% U - matrix with inputs, each column an input snapshot
% t - time vector
% rows - in how many rows should the figure be divided. The first row
%   always contains inputs, so the final rows will be nrows = rows + 1
% Sx - formatting for each input (e.g. 'c:', 'b-' etc.)
% Su - formatting for each state

%Checks:
% - if the desired number of rows and columns match the number of states
% - of X,U have same cols
% - if Sx and Su have enough elements

nx = size(X,1); %num of states
nu = size(U,1); %num of inputs
cols = ceil(nx/rows);
rows = rows + 1; % first row is the input plot

uhandles = zeros(cols, 1);
xhandles = zeros(nx, 1);

%Plot the input on the first row
for i = 1:cols
    subplot(rows, cols,i);
    hold on
    for j = 1:nu
        plot(t, U(j,:) ,Su{j});
        grid on
    end
    hold off
    uhandles(i) = gca;
end

%Plot the states
for i = 1:nx
    subplot(rows,cols,i+cols);
    plot(t,X(i,:), Sx{i});
    grid on
    xhandles(i) = gca;
end

