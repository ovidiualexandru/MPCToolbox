function h_disc = nonlinear_c2d(handle_nl_model, Ts, method)
%NONLINEAR_C2D Converts continuous-time nonlinear dynamic system to 
%discrete time nonlinear system.
%
%   H_DISC = NONLINEAR_C2D(HANDLE_NL_MODEL, Ts, METHOD) Gets the discrete
%   nonlinear model function handles H_DISC for the continous nonlinear
%   model described by the function handle HANDLE_NL_MODEL, using the
%   method from METHOD, assuming constant inputs over the time-step and
%   with sample time Ts.
%
%   Input arguments:
%   - handle_nl_model: the nonlinear model function handle. This
%   function must accept parameters compatible with the ode45 function
%   handle parameter (ODEFUN). This function must accept two arguments:
%   1. a time scalar value, which should not be used (for the euler
%   integration to work), in other words, it must be time invariant.
%   2. a (nx+nu)-by-1 vector in which the state in concatenated with the
%   input. The function should know to separate them and return the state
%   derivative in the first nx values of the return vector (which should be
%   (nx+nu)-by-1), padded with nu-by-1 zeros.
%   - Ts: the sampling time.
%   - method: the discretization method. Accepted values:
%   'ode45'     Integration using ode45. Is very accurate, but also very 
%               slow. For more information, type
%               help nonlinear_c2d>disc_ode
%   'euler'     Euler method for integration. A lot faster than ode45, but
%               the timestep needs to be small. For more information, type
%               help nonlinear_c2d>disc_euler
%
%   Output arguments:
%   - h_disc: the nonlinear discrete function handle.
    function x = disc_euler(x0,u)
%DISC_EULER Discrete Nonlinear Model with Euler discretization
%
%   X = DISC_EULER(X0,U) returns the state X starting from initial state X0
%   simulated over a timespte with constant input U.
%
%   Input arguments:
%   - X0: the starting (initial) state. Must be a nx-by-1 vector
%   - U: the input. Must be a nu-by-1 vector
%
%   Output arguments:
%   - X: the new state
        nx = length(x0);
        f = handle_nl_model([],[x0; u]);
        xd = f(1:nx);
        x = x0 + Ts*xd;
    end

    function x = disc_ode(x0,u)
%DISC_ODE Discrete Nonlinear Model with ODE45 integration.
%
%   X = DISC_ODE(X0, U) returns the state X starting from initial state X0, 
%   simulated over a timestep with constant input U.
%
%   Arguments:
%   - X0: the starting (initial) state. Must be a nx-by-1 vector
%   - U: the input. Must be a nu-by-1 vector
%
%   Output arguments:
%   - X: the new state
        nx = length(x0);
        [~, Yout] = ode45(handle_nl_model, [0 Ts], [x0; u]); %f(xk, uk)
        x = Yout(end, 1:nx)'; %get new real state
    end
switch method
    case 'euler'
        h_disc = @disc_euler;
    case 'ode45'
        h_disc = @disc_ode;
    otherwise
        error('Invalid method supplied. Please try ''euler'' or ''ode45''')
end
end
