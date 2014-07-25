function x = quanser_disc_nl2(x0,u,h)
%QUANSER_DISC_NL Discrete Nonlinear Model for the Quanser 3-DOF
%helicopter with ODE45 discretization.
%   x = QUANSER_DISC_NL(x0,u,h) returns the state x starting from
%   initial state x0, simulated over timestep h with constant input u.
%
%   Arguments:
%   - x0: the starting (initial) state. Must be a 6-by-1 vector
%   - u: the input. Must be a 2-by-1 vector
%   - h: the simulation time, or timestep. Must be a scalar.
%   Output arguments:
%   - x: the new state
%   Notes: must have quanser_cont_nl function in PATH.

[~, Yout] = ode45(@quanser_cont_nl2, [0 h], [x0; u]); %f(xk, uk)
x = Yout(end, 1:6)'; %get new real state
end

