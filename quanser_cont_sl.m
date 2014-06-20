function [A,B] = quanser_cont_sl(x0, u)
%% *Continous SL model for the Quanser helicopter*
% The model is taken from here: 
% <https://www.dropbox.com/s/lvvh5a2w9qkb2ll/chp_10.1007_978-94-007-6516-0_11.pdf>
% The state vector is defined ( <_d> meaning derived):
% x = [epsilon epsilon_d theta theta_d phi phi_d]';
% Parameters:
% - x0 : initial state.
% - u : a 2-by-1 vector containing the inputs
% - A,B - linearized model pair
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
