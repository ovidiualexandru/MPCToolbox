function [h_disc_nl_euler, h_disc_nl] = quanser_discrete(handle_cont_nl)
%QUANSER_DISCRETE Discretize nonlinear Quanser model
%   Explanation soon to come

    function x = disc_nl(x0,u,h)
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
        
        [~, Yout] = ode45(handle_cont_nl, [0 h], [x0; u]); %f(xk, uk)
        x = Yout(end, 1:6)'; %get new real state
    end

    function x = disc_nl_euler(x0,u,h)
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
        
        f = handle_cont_nl([],[x0; u]);
        xd = f(1:6);
        x = x0 + h*xd;
    end
h_disc_nl_euler = @disc_nl_euler;
h_disc_nl = @disc_nl;
end

