function [x_o, u_o] = affine_eq(A,B,g)
%AFFINE_EQ Summary of this function goes here
%   Detailed explanation goes here
nx = size(A,1);
u_o = linsolve([B], -g);
x_o = zeros(nx,1);
end

