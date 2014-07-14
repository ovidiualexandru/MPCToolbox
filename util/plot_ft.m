function plot_ft(FVAL, TEVAL, figtitle, fignumber)
%PLOT_FT Plot FVAL and TEVAL for simulation
%   PLOT_FT(FVAL, TEVAL, figtitle, fignumber) plots the values in FVAL
%   using TEVAL as time vector in figure <fignumber> and the title in
%   <figtitle>.
%
%   Arguments:
%   - FVAL: vector containing values
%   - TEVAL: time vector. Recommended: 1:length(FVAL)
%   - figtitle: the figure title
%   - fignumber: the figure number

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

