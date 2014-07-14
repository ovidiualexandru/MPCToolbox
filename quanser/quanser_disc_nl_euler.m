function x = quanser_disc_nl_euler(x0,u,h)
%QUANSER_DISC_NL_EULER Discrete Nonlinear Model for the Quanser 3-DOF
%helicopter with Euler discretization
%   x = QUANSER_DISC_NL_EULER(x0,u,h) returns the state x starting from
%   initial state x0, simulated over timestep h with constant input u.
%
%   Arguments:
%   - x0: the starting (initial) state. Must be a 6-by-1 vector
%   - u: the input. Must be a 2-by-1 vector
%   - h: the simulation time, or timestep. Must be a scalar.
%   Output arguments:
%   - x: the new state

f = quanser_cont_nl([],[x0; u]);
xd = f(1:6);
x = x0 + h*xd;
end

