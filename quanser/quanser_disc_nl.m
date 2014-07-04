function x = quanser_disc_nl(x0,u,h)
%QUANSER_DISC_NL Nonlinear Quanser model simulation
%   Uses ode45 to simulate over time h for fixed input u from initial state
%   x0. Uses the quanser_cont_nl function
[~, Yout] = ode45(@quanser_cont_nl, [0 h], [x0; u]); %f(xk, uk)
x = Yout(end, 1:6)'; %get new real state
end

