function [handle_nl_model, handle_sl_model] = quanser_model(parameters)
%QUANSER_MODEL Get the NL and SL continous model for the Quanser
%   This will be explained

%% Model data
% Set default values
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
niu_epsilon = 0.001; %kg*m^2/s
niu_theta = 0.001; %kg*m^2/s
niu_phi = 0.005; %kg*m^2/s
%Check if default values need to be changed
if isfield(parameters,'Jtheta') && ~isempty(parameters.Jtheta)
    Jtheta = parameters.Jtheta;
end
if isfield(parameters,'Jphi') && ~isempty(parameters.Jphi)
    Jphi = parameters.Jphi;
end
if isfield(parameters,'La') && ~isempty(parameters.La)
    La = parameters.La;
end
if isfield(parameters,'Lc') && ~isempty(parameters.Lc)
    Lc = parameters.Lc;
end
if isfield(parameters,'Ld') && ~isempty(parameters.Ld)
    Ld = parameters.Ld;
end
if isfield(parameters,'Le') && ~isempty(parameters.Le)
    Le = parameters.Le;
end
if isfield(parameters,'Lh') && ~isempty(parameters.Lh)
    Lh = parameters.Lh;
end
if isfield(parameters,'Mf') && ~isempty(parameters.Mf)
    Mf = parameters.Mf;
end
if isfield(parameters,'Mb') && ~isempty(parameters.Mb)
    Mb = parameters.Mb;
end
if isfield(parameters,'Mc') && ~isempty(parameters.Mc)
    Mc = parameters.Mc;
end
if isfield(parameters,'Km') && ~isempty(parameters.Km)
    Km = parameters.Km;
end
if isfield(parameters,'niu_epsilon') && ~isempty(parameters.niu_epsilon)
    niu_epsilon = parameters.niu_epsilon;
end
if isfield(parameters,'niu_theta') && ~isempty(parameters.niu_theta)
    niu_theta = parameters.niu_theta;
end
if isfield(parameters,'niu_phi') && ~isempty(parameters.niu_phi)
    niu_phi = parameters.niu_phi;
end
%Fixed values
g = 9.81; %m/s^2;
deltaa = atan((Ld+Le)/La);
deltac = atan(Ld/Lc);
deltah = atan(Le/Lh);
%% Compute coefficients
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
%% Nonlinear function
    function f = nl(t, y)
        %QUANSER_CONT_NL Continous nonlinear model for the Quanser 3-DOF
        %helicopter.
        %   f = QUANSER_CONT_NL(t, y) compute the states derivative, f using the
        %   initial state x0 and input u concatenated into vector y = [x0; u]
        %
        %   Arguments:
        %   - t : time-instant (required by ode45), not used in function.
        %   - y : an 8-by-1 vector with the initial state and inputs. y = [x0; u]
        %   Output arguments:
        %   - f : the derived state, padded with zeros (assuming constant inputs).
        %       f = [xd; zeros(2,1)]; where xd = F + G*u
        %
        %   The model is taken from:
        %https://www.dropbox.com/s/lvvh5a2w9qkb2ll/chp_10.1007_978-94-007-6516-0_11.pdf
        %   The state vector is defined ( <_d> meaning derived):
        %              x = [epsilon epsilon_d theta theta_d phi phi_d]';
        
        x0 = y(1:6);
        u = y(7:8);
        %% Model
        epsilon = x0(1);
        epsilon_d = x0(2);
        theta = x0(3);
        theta_d = x0(4);
        phi = x0(5);
        phi_d = x0(6);
        %% Model matrices
        f1 = epsilon_d;
        f2 = p(1)*cos(epsilon*pi/180) + p(2)*sin(epsilon*pi/180) + p(3)*epsilon_d;
        f3 = theta_d;
        f4 = p(5)*cos(theta*pi/180) + p(6)*sin(theta*pi/180) + p(7)*theta_d;
        f5 = phi_d;
        f6 = p(9)*phi_d;
        F = [ f1; f2; f3; f4; f5; f6];
        
        g2 = p(4)*cos(theta*pi/180);
        g4 = p(8);
        g6 = p(10)*sin(theta*pi/180);
        G1 = [0; g2; 0; g4; 0; g6];
        G2 = [0; g2; 0; -g4; 0; g6];
        G = [G1 G2];
        %% Calculate derived state
        xd = F + G*u;
        f = [xd; zeros(2,1)];
    end
%% SL function
    function [A,B,g] = sl(x0, u)
        %SL_MODEL Continous Successive liniarization model for the Quanser
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
        f = nl([], [x0; u]);
        xd = f(1:6);
        g = xd - A*x0 - B*u;
    end
%% Return handles
handle_nl_model = @nl;
handle_sl_model = @sl;
end