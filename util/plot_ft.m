function plot_ft(FVAL, TEVAL, figtitle, fignumber)
%PLOT_FT Plot FVAL and TEVAL for simulation
%
%   PLOT_FT(FVAL, TEVAL, FIGTITLE, FIGNUMBER) plots the values in FVAL 
%   and TEVAL in figure FIGNUMBER and the title in FIGTITLE.
%
%   Input arguments:
%   - FVAL: vector containing values, usually from the objective function.
%   - TEVAL: vector containing the evaluation time for each value in FVAL.
%   - FIGTITLE: the figure title.
%   - FIGNUMBER: the figure number.

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

