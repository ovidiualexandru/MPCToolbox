function quanser_phase_plot(X,figtitle)
%% Quanser phase plot
% Call this for a nice phase plot for each state pair.
% Examples: quanser_mpc, quanser_lqr
%% Configuration
titles = {'Elevation $\epsilon$'; 'Pitch $\theta$'; 'Travel $\phi$'};
%% Figure initialization
rows = 1;
cols = 3;
figure(2);
clf;
whitebg([1 1 1]);
figure(2);
set(gcf, 'Name',figtitle);
%% Phase-plot
for i = 1:3
    k = 2*i - 1; %state index
    subplot(rows, cols, i);
    plot(X(k,:), X(k+1,:), 'b-');
    xlabel('Angle [deg]');
    ylabel('Angular speed [deg/s]');
    title(titles{i},'Interpreter','latex');
    grid on
end
