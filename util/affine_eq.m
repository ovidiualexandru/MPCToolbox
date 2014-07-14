function [x_o, u_o] = affine_eq(A,B,g)
%AFFINE_EQ Get the state-input pair that eliminates the the affine term.
%   [x_o, u_o] = affine_eq(A,B,g). Find x_o and u_o so that the the affine
%   term in the system dynamic is equal to 0.
%   Given the system dynamic equation:
%                   x_dot = A*x + B*u + g
%   Find a pair (x_o, u_o) so that xbar = x - x_o and ubar = u - u_o gives:
%                   xbar_dot = A*xbar + B*ubar.
%   Note: in this implementation, x_o is always 0.
%
%   Arguments:
%   - A, B: the state-space matrices describing the system dynamic.
%   - g: the affine term in the system dynamic.
%   Output arguments:
%   - x_o, u_o: the state-input pair that eliminates the affine term, they
%   must be added to x and u respectively to get the true state and input.
nx = size(A,1);
u_o = linsolve([B], -g);
x_o = zeros(nx,1);
end

