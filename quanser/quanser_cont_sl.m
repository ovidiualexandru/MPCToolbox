function [A,B,g] = quanser_cont_sl(x0, u)
%QUANSER_CONT_SL Continous Successive liniarization model for the Quanser
%3-DOF model.
%   [A,B,g] = QUANSER_CONT_SL(x0, u) gets the liniarized model pair (A,B,g)
%   from initial state x0 and input u.
%
%   Arguments:
%   - x0 : initial state, a 6-by-1 vector.
%   - u : a 2-by-1 vector containing the inputs
%   Output arguments:
%   - A,B : linearized model pair
%   - g : the affine term in the linear system dynamic:
%                       x_dot = A*x + B*u + g, where
%                      g = f(x0,u0) - A*x0 - B*u0, and
%      f(x0,u0) - x_dot from the NL contionous model (quanser_cont_nl)
%
%   The model is taken from: 
%https://www.dropbox.com/s/lvvh5a2w9qkb2ll/chp_10.1007_978-94-007-6516-0_11.pdf
%   The state vector is defined ( <_d> meaning derived):
%              x = [epsilon epsilon_d theta theta_d phi phi_d]';
%   Notes: must have quanser_cont_nl function in PATH.

%% Model params
epsilon = x0(1);
epsilon_d = x0(2);
theta = x0(3);
theta_d = x0(4);
phi = x0(5);
phi_d = x0(6);
Vf = u(1);
Vb = u(2);
p = quanser_params();
%% Model matrices
a21 = -p(1) * sin(epsilon*pi/180) + p(2) * cos(epsilon*pi/180);
a22 = p(3);
a23 = -p(4) * sin(theta*pi/180) * (Vf + Vb);
a43 = -p(5) * sin(theta*pi/180) + p(6) * cos(theta*pi/180);
a44 = p(7);
a63 = p(10) * cos(theta*pi/180) * (Vf+Vb);
a66 = p(9);
A = [ 0  , 1  , 0  , 0  , 0  , 0  ;
      a21, a22, a23, 0  , 0  , 0  ;
      0  , 0  , 0  , 1  , 0  , 0  ;
      0  , 0  , a43, a44, 0  , 0  ;
      0  , 0  , 0  , 0  , 0  , 1  ;
      0  , 0  , a63, 0  , 0  , a66;
    ];

g2 = p(4)*cos(theta*pi/180);
g4 = p(8);
g6 = p(10)*sin(theta*pi/180);
G1 = [0; g2; 0; g4; 0; g6];
G2 = [0; g2; 0; -g4; 0; g6];
B = [G1 G2];
f = quanser_cont_nl([], [x0; u]);
xd = f(1:6);
g = xd - A*x0 - B*u;
