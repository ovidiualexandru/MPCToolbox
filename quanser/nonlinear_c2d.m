function [h_disc_euler, h_disc_ode] = nonlinear_c2d(handle_nl_model, Ts)
%NONLINEAR_C2D Discretize nonlinear model
%   [handle_nl_model, handle_sl_model] = NONLINEAR_C2D(handle_nl_model, Ts) 
%   Gets the discretized nonlinear model as function handles, using euler 
%   or ode45 integration methods, assuming constant inputs over the
%   time-step, with sample time Ts.
%
%   Input arguments:
%       - handle_nl_model: the nonlinear model function handle. This
%           function must accept parameters compatible with the ode45
%           function handle parameter (ODEFUN). This function must accept 
%           two arguments:
%               1. a time scalar value, which should not be used (for the
%                   euler integration to work).
%               2. a (nx+nu)-by-1 vector in which the state in concatenated
%                   with the input. The function should know to separate
%                   them and return the state derivative in the first nx
%                   values of the return vector (which should be
%                   (nx+nu)-by-1), padded with nu-by-1 zeros.
%       - Ts: the sampling time.
%   Output arguments:
%       - h_disc_nl_euler: the nonlinear discrete function handle, using
%           euler method integration. For more information, type:
%           help nonlinear_c2d>disc_euler
%       - h_disc_ode: the nonlinear discrete function handle, using
%           ode45 integration. For more information, type:
%           help nonlinear_c2d>disc_ode
    function x = disc_euler(x0,u)
        %DISC_EULER Discrete Nonlinear Model with Euler discretization
        %   x = DISC_EULER(x0,u,h) returns the state x starting from
        %   initial state x0, simulated over timestep h with constant input u.
        %
        %   Arguments:
        %   - x0: the starting (initial) state. Must be a nx-by-1 vector
        %   - u: the input. Must be a nu-by-1 vector
        %   Output arguments:
        %   - x: the new state
        nx = length(x0);
        f = handle_nl_model([],[x0; u]);
        xd = f(1:nx);
        x = x0 + Ts*xd;
    end

    function x = disc_ode(x0,u)
        %DISC_ODE Discrete Nonlinear Model with ODE45 integration.
        %   x = DISC_ODE(x0,u,h) returns the state x starting from
        %   initial state x0, simulated over timestep h with constant input
        %   u.
        %
        %   Arguments:
        %   - x0: the starting (initial) state. Must be a nx-by-1 vector
        %   - u: the input. Must be a nu-by-1 vector
        %   Output arguments:
        %   - x: the new state
        nx = length(x0);
        [~, Yout] = ode45(handle_nl_model, [0 Ts], [x0; u]); %f(xk, uk)
        x = Yout(end, 1:nx)'; %get new real state
    end
h_disc_euler = @disc_euler;
h_disc_ode = @disc_ode;
end
