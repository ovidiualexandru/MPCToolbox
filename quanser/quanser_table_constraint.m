function c = quanser_table_constraint(x)
%QUANSER_TABLE_CONSTRAINT Nonlinear constraint on the elevation and pitch
%angles so that the helicopter does not touch the table.
%
%   C = QUANSER_TABLE_CONSTRAINT(X) returns a value in C which contains
%   information about the lowest point of the helicopter relative to the
%   table. This value should always be <= 0. This function is implemented 
%   to be compatible with the nonlinear state constraint (fcu) function in 
%   nmpc_fullspace.
%
%   Input arguments:
%   - X: a 6-by-1 vector with the current state of the helicopter.
%
%   Output arguments:
%   - C: a 6-by-1 vector with the first value representing the constraint, 
%   padded with zeros.
%
%   Notes: The distances used are default values from documentation. If the
%   model is changed this function will not work properly.

%% Variables
epsilon = x(1); %elevation angle
theta = x(3); %pitch angle
%% Fixed values
% TODO : check the h and l variables
La = 0.62; %[m] Distance from the pivot point to the helicopter body
Lh = 0.177; %[m] Distance from the pitch axis to either motor
l = 0.1;% [m] Distance from motor to edge of helicopter
h = 0.3;%[m] Height of the pivot point
w = Lh + l; %[m] Total distance from pitch axis to edge of helicopter
%% Computations
c = zeros(6,1);
c(1) = La*sin(-epsilon*pi/180) + w*abs(sin(theta*pi/180)) - h;
end

