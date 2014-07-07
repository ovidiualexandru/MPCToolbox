function [x_o, u_o] = affine_eq(A,B,g)
%AFFINE_EQ Summary of this function goes here
%   Detailed explanation goes here
nx = size(A,1);
n = linsolve([A B], -g);
x_o = n(1:nx);
u_o = n(nx + 1:end);
end

