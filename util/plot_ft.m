function plot_ft(FVAL, TEVAL, figtitle, fignumber)
%PLOT_FT Plot FVAL and TEVAL for simulation
%   Make a nice plot for FVAL (Cost function from solver) and TEVAL (time 
%   needed to compute a command).

N = length(FVAL);
t = 1:N;
figure(fignumber);
clf;
set(gcf, 'Name',figtitle);
whitebg([0 0 0]);

subplot(2, 1, 1);
plot(t, FVAL);
title('Solver cost values');
xlabel('[k]');
ylabel('FVAL');
grid on

subplot(2, 1, 2);
stem(t,TEVAL.*1000);
title('Solver time');
xlabel('[k]');
ylabel('Eval time [ms]');
grid on
end

