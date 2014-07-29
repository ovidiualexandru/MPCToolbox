function p = quanser_params(param)
%QUANSER_PARAMS returns the Quanser 3-DOF helicopter nonlinear model
%parameters.
%
%   P = QUANSER_PARAMS() returns in vector P the 10 coefficients for the 
%   model using default values.
%
%   P = QUANSER_PARAMS(PARAM) returns in vector P the 10 coefficients for
%   the model using the fields in structure PARAM.
%
%   Input arguments:
%   - PARAM: structure containing the model parameters as fields. These 
%   parameters are enumerated below. If the structure does not have one of
%   these fields, the default value for that field is used.
%   PARAM.Jepsilon: Moment of inertia about the elevation axis.
%       Default value: 0.86 kg*m^2
%   PARAM.Jtheta: Moment of inertia about the pitch axis. 
%       Default value: 0.044 kg*m^2
%   PARAM.Jphi: Moment of inertia about the travel axis.
%       Default value: 0.82 kg*m^2
%   PARAM.La: Distance from the pivot point to the helicopter body.
%       Default value: 0.62 m
%   PARAM.Lc: Distance from the pivot point to the counter-weight.
%       Default value: 0.44 m
%   PARAM.Ld: Length of pendulum for the elevation axis.
%       Default value: 0.05 m
%   PARAM.Le: Length of pendulum for the pitch axis.
%       Default value: 0.02 m
%   PARAM.Lh: Distance from the pitch axis to either motor.
%       Default value: 0.177 m
%   PARAM.Mf: Mass of the front section of the helicopter.
%       Default value: 0.69 kg
%   PARAM.Mb: Mass of the rear section.
%       Default value: 0.69 kg
%   PARAM.Mc: Mass of the counter-weight.
%       Default value: 1.69 kg
%   PARAM.Km: The voltage-force gain.
%       Default value: 0.5 N/V
%   PARAM.niu_epsilon: Coefficient of viscous friction - elevation.
%       Default value: 0.001 kg*m^2/s
%   PARAM.niu_theta: Coefficient of viscous friction - pitch.
%       Default value: 0.001 kg*m^2/s
%   PARAM.niu_phi: Coefficient of viscous friction - travel.
%       Default value: 0.005 kg*m^2/s
%
%   Output arguments:
%   - P: a 10-by-1 vector containing the p1..p10 coefficients, used in 
%   pages 3 and 5 in
%https://www.dropbox.com/s/lvvh5a2w9qkb2ll/chp_10.1007_978-94-007-6516-0_11.pdf

%% Argument processing
if nargin == 0
    param = [];
end
%% Model data
Jepsilon = 0.86; %[kg*m^2] Moment of inertia about the elevation axis
Jtheta = 0.044; %[kg*m^2] Moment of inertia about the pitch axis
Jphi = 0.82; %%[kg*m^2] Moment of inertia about the travel axis
La = 0.62; %[m] Distance from the pivot point to the helicopter body
Lc = 0.44; %[m] Distance from the pivot point to the counter-weight
Ld = 0.05; %[m] Length of pendulum for the elevation axis
Le = 0.02; %[m] Length of pendulum for the pitch axis
Lh = 0.177; %[m] Distance from the pitch axis to either motor
Mf = 0.69; %[kg] Mass of the front section of the helicopter
Mb = 0.69; %[kg] Mass of the rear section
Mc = 1.69; %[kg] Mass of the counter-weight
Km = 0.5; %[N/V]
niu_epsilon = 0.001;%[kg*m^2/s] Coefficient of viscous friction - elevation
niu_theta = 0.001; %[kg*m^2/s] Coefficient of viscous friction - pitch
niu_phi = 0.005; %[kg*m^2/s] Coefficient of viscous friction - travel
%Check if default values need to be changed
if isfield(param,'Jtheta') && ~isempty(param.Jtheta)
    Jtheta = param.Jtheta;
end
if isfield(param,'Jphi') && ~isempty(param.Jphi)
    Jphi = param.Jphi;
end
if isfield(param,'La') && ~isempty(param.La)
    La = param.La;
end
if isfield(param,'Lc') && ~isempty(param.Lc)
    Lc = param.Lc;
end
if isfield(param,'Ld') && ~isempty(param.Ld)
    Ld = param.Ld;
end
if isfield(param,'Le') && ~isempty(param.Le)
    Le = param.Le;
end
if isfield(param,'Lh') && ~isempty(param.Lh)
    Lh = param.Lh;
end
if isfield(param,'Mf') && ~isempty(param.Mf)
    Mf = param.Mf;
end
if isfield(param,'Mb') && ~isempty(param.Mb)
    Mb = param.Mb;
end
if isfield(param,'Mc') && ~isempty(param.Mc)
    Mc = param.Mc;
end
if isfield(param,'Km') && ~isempty(param.Km)
    Km = param.Km;
end
if isfield(param,'niu_epsilon') && ~isempty(param.niu_epsilon)
    niu_epsilon = param.niu_epsilon;
end
if isfield(param,'niu_theta') && ~isempty(param.niu_theta)
    niu_theta = param.niu_theta;
end
if isfield(param,'niu_phi') && ~isempty(param.niu_phi)
    niu_phi = param.niu_phi;
end
%Fixed values
g = 9.81; %m/s^2;
deltaa = atan((Ld+Le)/La);
deltac = atan(Ld/Lc);
deltah = atan(Le/Lh);
%% Return coefficients
p = zeros(10,1);
p(1) = ( -(Mf + Mb) * g * La + Mc * g * Lc) / Jepsilon;
p(2) = ( -(Mf + Mb) * g * La * tan(deltaa) + Mc * g * Lc * tan(deltac) )...
    / Jepsilon;
p(3) = -niu_epsilon / Jepsilon;
p(4) = Km * La / Jepsilon;
p(5) = (-Mf + Mb) * g * Lh / Jtheta;
p(6) = -(Mf + Mb) * g * Lh * tan(deltah) / Jtheta;
p(7) = -niu_theta / Jtheta;
p(8) = Km * Lh / Jtheta;
p(9) = -niu_phi / Jphi;
p(10) = -Km * La / Jphi;
