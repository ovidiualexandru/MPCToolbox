clear
addpath('./quanser');
addpath('./util');
%% System initialization
x0 = [0; 0; 0; 0; 0; 0]; %Initial state
u0 = [1.8; 1.8]; % [Vf Vb] initial inputs
h = 0.1; % s - sampling time
nu = 2;
nx = 6;
L = 3; % Simulation progress update rate
Nc = 3; % Control and prediction horizon
N = 600; % Simulation size
%% Cost matrices and constraints
Q = diag([5, .001, .1, .001, 0, 0],0);
R = diag([.01, .01],0);
%state constraints, positive and negative
dx = [ 30,  100,  100,  50,  inf,  inf;
      -30, -100, -100, -50, -inf, -inf];
%input constraints
du = [ 5,  5;
       0,  0];
 %% Model generation
% Set model coefficients. Leave empty for default value
mpc_param= []; % Use nominal model
%Get MPC continous model
mpc_sl = quanser_model('sl', mpc_param);
%% Problem initialization
% Fields not defined here are not fixed and are defined below
problem = struct;
problem.Q = Q;
problem.R = R;
problem.Nc = Nc;
%% Solver initialization
X = zeros(nx, N); %save all states, for plotting
U = zeros(nu, N); %save all inputs
XREF = zeros(nx, N); %save reference
UREF = zeros(nu, N); %save input reference
FVAL = zeros(1, N); %save cost value
TEVAL = zeros(1, N); %save calculation time
x = x0;
xr = x0; % 'real' x
u = u0;
uref = u0;
%% Setup
statedef;
% URI used to connect to the server model.
uri = 'shmem://foobar:1';
% Use blocking I/O. Do not change this value.
nonblocking = false;
% Connecting to the server using the specified URI
fprintf(1, 'Connecting to the server at URI: %s\n', uri);
stream = stream_connect(uri, nonblocking);
%Messages
fprintf(1, 'Connected to server.\n\n');
%% Process control loop
ue = []; %input estimated solution
try
    for i=1:N
    %% Update SL Model
        tic;
        %% Receive data
        value = stream_receive_double_array(stream,8);
        x = zeros(6,1);
        xref = zeros(6,1);
        x(1) = value(3); %Elevation
        x(2) = value(4); %Elevation Rate
        x(3) = value(5); %Pitch
        x(4) = value(6); %Pitch Rate
        x(5) = value(7); %Travel
        x(6) = value(8); %Travel Rate
        xref(3) = value(1); %Pitch Ref
        xref(1) = value(2); %Elevation Ref
        %Do work
        if mod(i,L) == 0 || i == 1
            [A,B,g] = mpc_sl(x,u); %recalculate (A,B,g)
            [x_o, u_o] = affine_eq(A,B,g);
            du_bar = du - repmat(u_o',2,1);
            dx_bar = dx - repmat(x_o',2,1);
            Ad = eye(nx) + h*A;
            Bd = h*B;
            problem.A = Ad;
            problem.B = Bd;
            problem.du = du_bar;
            problem.dx = dx_bar;
            fprintf('%d ', i);
            if mod(i,20*L) == 0
                fprintf('\n');
            end
        end
    %% Get next command
    xbar = x - x_o;
    idif = Nc - 1;
    if i + Nc > N
        idif = N - i;
    end
    urefbar = uref - u_o;
    xrefbar = xref - x_o;
    problem.uref = urefbar;
    problem.xref = xrefbar;
    problem.uprev = ue;
    problem.x0 = xbar;
    [ue, Xe,fval,EXITFLAG, OUTPUT] = lmpc_condensed(problem);
    if EXITFLAG < 0
        fprintf('Iteration %d\n',i)
        stream_close(stream);
        fprintf(1, 'Connection closed\n');
        error('Solver error \n');
    end
    ubar = ue(:,1); %use only the first command in the sequence
    u = ubar + u_o;
    teval = toc;
    %% Data logging
    X(:,i) = x; % save states
    U(:,i) = u; % save inputs
    XREF(:,i) = xref;
    FVAL(i) = fval;
    TEVAL(i) = teval;
    %% Send data
    % Store a double value in the stream send buffer
    us = struct;
    us.Vf = u(1);
    us.Vb = u(2);
    stream_send_array(stream, us);
    % Flush the send buffer to the underlying communications channel
    stream_flush(stream);
    if isempty(value) % then the server closed the connection gracefully
        fprintf(1, '\nServer has closed the connection.\n');
        break;
    end
    end
    % Once the Esc key is pressed, close the stream handle used for
    % communications.
    fprintf(1, '\nShutting down the client...\n');
    stream_close(stream);
    fprintf(1, 'Connection closed\n');
catch
    err = lasterror;
    fprintf(1, '\n%s.\nShutting down the client from catch...\n', err.message);
    stream_close(stream);
    fprintf(1, 'Connection closed\n');
    rethrow(err);
end
%% Plotting
quanser_plot(X, U, dx, du,['titlu' ' Quanser Plot'], 7, XREF);
quanser_phase_plot(X, ['titlu' ' Quanser Phase-Plot'], 8, XREF);
plot_ft(FVAL, TEVAL, ['titlu' ' Quanser Performance'], 9);