function [x_o, u_o] = affine_eq(A,B,g)
%AFFINE_EQ Get the state-input pair that eliminates the the affine term in
%a linear system dynamic
%
%   [X_O, U_O] = AFFINE_EQ(A,B,G). Find X_O and U_O so that the the affine
%   term in the system dynamic is equal to 0.
%   Given the system dynamic equation:
%                   x_dot = A*x + B*u + g
%   Find a pair (x_o, u_o) so that xbar = x - x_o and ubar = u - u_o gives:
%                   xbar_dot = A*xbar + B*ubar.
%   Note: in this implementation, x_o is always 0.
%
%   Input arguments:
%   - A, B: the state-space matrices describing the system dynamic.
%   - G: the affine term in the system dynamic.
%
%   Output arguments:
%   - X_O, U_O: the state-input pair that eliminates the affine term, they
%   must be added to x and u respectively to get the true state and input.
nx = size(A,1);
u_o = linsolve([B], -g);
x_o = zeros(nx,1);
end

