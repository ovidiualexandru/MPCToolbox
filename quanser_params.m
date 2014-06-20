function p = quanser_params()
%% *Get Quanser model parameters*
% The model is taken from here: 
% <https://www.dropbox.com/s/lvvh5a2w9qkb2ll/chp_10.1007_978-94-007-6516-0_11.pdf>
% Output:
% - p: vector containing the p1...p10 coefficients, used in page 3 and 5.
%% Model data
Jepsilon = 0.86; %kg*m^2
Jtheta = 0.044; %kg*m^2
Jphi = 0.82; %kg*m^2
La = 0.62; %m
Lc = 0.44; %m
Ld = 0.05; %m
Le = 0.02; %m
Lh = 0.177; %m
Mf = 0.69; %kg
Mb = 0.69; %kg
Mc = 1.69; %kg
Km = 0.5; %N/V
g = 9.81; %m/s^2;
niu_epsilon = 0.001; %kg*m^2/s
niu_theta = 0.001; %kg*m^2/s
niu_phi = 0.005; %kg*m^2/s
deltaa = atan((Ld+Le)/La);
deltac = atan(Ld/Lc);
deltah = atan(Le/Lh);
%% Return coefficients
p = zeros(10,1);
p(1) = ( -(Mf + Mb) * g * La + Mc * g * Lc) / Jepsilon;
p(2) = ( -(Mf + Mb) * g * La * tan(deltaa) + Mc * g * Lc * tan(deltac) ) / Jepsilon;
p(3) = -niu_epsilon / Jepsilon;
p(4) = Km * La / Jepsilon;
p(5) = (-Mf + Mb) * g * Lh / Jtheta;
p(6) = -(Mf + Mb) * g * Lh * tan(deltah) / Jtheta;
p(7) = -niu_theta / Jtheta;
p(8) = Km * Lh / Jtheta;
p(9) = -niu_phi / Jphi;
p(10) = -Km * La / Jphi;