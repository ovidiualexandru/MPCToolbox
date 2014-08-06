function [u, X, FVAL, EXITFLAG, OUTPUT] = nmpc_fullspace(problem)
%NMPC_FULLSPACE Compute the input sequence and predicted output using
%Nonlinear MPC fullspace (simultaneous) approach.
%
%   [U, X, FVAL, EXITFLAG, OUTPUT] = NMPC_FULLSPACE(PROBLEM). Compute the
%   input sequence using a simultaneous Non-Linear MPC formulation over Nc
%   samples for the nonlinear discrete system fd(x,u) with (Q,R) as weight
%   matrices, starting from state x0 and constraints on inputs desribed by
%   du and on states by dx. The fields of PROBLEM are described below. Some
%   fields are mandatory, others are optional. The optional fields can be
%   omitted or set to empty array.
%   Mathematical formulation of the problem:
%
%              Nc-1
%              --
%              \
%          min /   (x - xref )'Q(x - xref ) + (u - uref )'R(u - uref )
%           u  --    k      k     k      k      k      k     k      k
%              k=0
%                         s.t.:
%                             x(k+1) = fd(x(k), u(k))
%                             ul < u(k) < uu (input bounds)
%                             xl < x(k) < xu (state bounds)
%                             Gu*u(k) <= pu
%                             Gx*x(k) <= px, and
%                             du = [ul'; uu']
%                             dx = [xl'; xu']
%
%   Mandatory fields:
%   - fd: the discrete nonlinear model function handle
%   - Q,R: the weighting matrices in the cost function.
%   - Nc: the control horizon.
%   - x0: the current( initial) state of the system.
%   - du, dx: the constraint vectors for inputs and states. du is a 2-by-nu
%   matrix containing constraints for inputs. First line is upper bound,
%   second is lower bound for each input. dx is a 2-by-nx matrix with
%   constraints for states. If the input/state has no lower bound, set it's
%   corresponding value to -Inf. Conversely, if the input/state has no
%   upper bound, set to Inf.
%       nu - number of inputs, nx - number of states.
%
%   Optional fields:
%   - xref: the desired( reference) state. Must have nx lines, but can have
%   number of columns in the range [1, Nc].
%   - uref: the reference input (stabilizing input). Must have nu lines,
%   but can have number of columns in the range [1, Nc].
%   - xprev: the previously obtained state prediction. This will be used as
%   a starting point for the algorithm, along with uprev.
%   - uprev: the previously obtained input solution. This will be used to 
%   'warm start' the algorithm. This should be the u obtained at a previous 
%   step.
%   - Gx, px: constraint matrix and vector for states. Gx is a nrx-by-nx
%   matrix of combinations of constraints on the inputs and pu is a nrx
%   vector.
%       nrx - number of constraints on the combinations of states
%       Gx*x <= px
%   - Gu, pu: constraint matrix and vector for inputs. Gu is a nru-by-nu
%   matrix of combinations of constraints on the inputs and pu is a nru
%   vector.
%       nru - number of constraints on the combinations of inputs
%       Gu*u <= pu
%   - Gpec, ppec: constraint matrix and vector for the Persistent Exciting
%   Condition (PEC). These are constraints for the inputs that apply only
%   for the first prediction in the horizon. This does not guarantee the
%   PEC, but help can be used to implement it by applying the constraints
%   on the first input in the sequence. In all other aspects, is the same
%   as Gu and pu.
%   - fcu: nonlinear constraints on inputs. Is a function handle that
%   takes an input vector and returns a nu-by-1 vector with the constraint
%   values. Calling: cu = fcu(u). The constraints should be formulated so
%   that cu <= 0. Restrictions: only nu nonlinear constraints on inputs.
%   - fcx: nonlinear constraints on states. Is a function handle that
%   takes a state vector and returns a nx-by-1 vector with the constraint
%   values. Calling: cx = fcx(x). The constraints should be formulated so
%   that cx <= 0. Restrictions: only nx nonlinear constraints on inputs.
%
%   Output arguments:
%   - U: a nu-by-Nc matrix of computed inputs. U(:,1) must be used.
%   - X: a nx-by-Nc matrix of predicted states.
%   - FVAL: the object function value given by the numerical solver, 
%       fmincon.
%   - EXITFLAG: the exitflag from the solver. See 'help fmincon' for 
%       details. EXITFLAG is > 0 if a solution has been found.
%   - OUTPUT: the output from the solver. See 'help fmincon' for details.

%% Extract parameters from struct
% Required fields
handle_nlmodeld = problem.fd;
Q = problem.Q;
R = problem.R;
Nc = problem.Nc;
du = problem.du;
dx = problem.dx;
x0 = problem.x0;
% Optional fields: references and previous solutions
if isfield(problem, 'xref')
    xref = problem.xref;
else
    xref = [];
end
if isfield(problem,'uref')
    uref = problem.uref;
else
    uref = [];
end
if isfield(problem, 'xprev')
    Xprev = problem.xprev;
else
    Xprev = [];
end
if isfield(problem, 'uprev')
    uprev = problem.uprev;
else
    uprev = [];
end
% Optional fields: inputs/states combination constraints
if isfield(problem, 'Gx')
    Gx = problem.Gx;
else
    Gx = [];
end
if isfield(problem, 'px')
    px = problem.px;
else
    px = [];
end
if isfield(problem, 'Gu')
    Gu = problem.Gu;
else
    Gu = [];
end
if isfield(problem, 'pu')
    pu = problem.pu;
else
    pu = [];
end
if isfield(problem, 'Gpec')
    Gpec = problem.Gpec;
else
    Gpec = [];
end
if isfield(problem, 'ppec')
    ppec = problem.ppec;
else
    ppec = [];
end
if isfield(problem, 'fcx')
    fcx = problem.fcx;
else
    fcx = [];
end
if isfield(problem, 'fcu')
    fcu = problem.fcu;
else
    fcu = [];
end
%% Argument processing
nu = size(du,2); %number of inputs
nx = size(dx,2); %number of states
if isempty(xref)
    xref = zeros(nx, Nc);
end
if isempty(uref)
    uref = zeros(nu, Nc);
end
difx = Nc - size(xref,2);
difu = Nc - size(uref, 2);
% If xref does not have enough columns, append the last column difx times
if difx > 0
    xref = [xref, repmat(xref(:,end), [1 difx])];
end
% For uref same as for xref above
if difu > 0
    uref = [uref, repmat(uref(:,end), [1 difu])];
end
if isempty(Xprev)
    %if there is no previous solution, use the reference as a start point
    Xprev = xref;
else
    Xprev = [Xprev(:,1:end-1), xref(:,end)]; %shift the previous solution
end
if isempty(uprev)
    %if there is no previous solution, use the reference as a start point
    uprev = uref;
else
    uprev = [uprev(:,1:end-1), uref(:,end)]; %shift the previous solution
end
Zprev = [uprev; Xprev];
if ~isa(handle_nlmodeld, 'function_handle')
    error('fd must be a function handle.');
end
if ~isempty(fcx)
    if ~isa(fcx, 'function_handle')
        error('fcx must be a function handle.');
    end
end
if ~isempty(fcu)
    if ~isa(fcu, 'function_handle')
        error('fcu must be a function handle.');
    end
end
%% Nonlinear constraints function
    function [C,Ceq] = nonlconfunc(z)
        C = zeros(size(z));
        Xl = reshape(z, nu+nx,[]);
        ul = Xl(1:nu,:);
        Xl = Xl(nu+1:end,:);
        x = x0;
        Ceq = zeros(size(z));
        for i = 1:Nc
            b = Xl(:,i);
            x = handle_nlmodeld(x, ul(:,i));
            dif = b - x; %x_{k+1} - f(x_k)
            ceq = [zeros(nu, 1); dif]; %inputs have no eq constraints
            idx_start = (i-1)*(nx+nu)+1;
            idx_end = i*(nx+nu);
            Ceq(idx_start:idx_end) = ceq;
            if  ~isempty(fcu)
                C(idx_start:idx_start+nu) = fcu(ul(:,i));
            end
            if ~isempty(fcx)
                C(idx_start+nu+1:idx_start+nu+nx) = fcx(Xl(:,i));
            end
        end
    end
%% QP definition
ubx = dx(1,:)';
lbx = dx(2,:)';
ubu = du(1,:)';
lbu = du(2,:)';
lb = [lbu; lbx];
ub = [ubu; ubx];
LB = repmat(lb, Nc, 1);
UB = repmat(ub, Nc, 1);
Qsmall = blkdiag(R,Q);
Q_hat = Qsmall;
for i = 1:Nc-1
    %Add another element to the block diagonal matrices
    Q_hat = blkdiag(Q_hat, Qsmall);
end
zsmall = [ uref; xref];
zref = zsmall(:);
q = -Q_hat*zref;
z0 = Zprev(:);
%% Gbar, pbar
if isempty(Gu)
    Gu = zeros(1, nu);
    pu = 0;
end
if isempty(Gx)
    Gx = zeros(1, nx);
    px = 0;
end
Gbar = [];
Gbarsmall = blkdiag(Gu, Gx);
for i = 1:Nc
    Gbar = blkdiag(Gbar, Gbarsmall);
end
pbarsmall = [pu; px];
pbar = repmat(pbarsmall, [Nc 1]);
if ~isempty(Gpec)
    % Add the pec constraint
    nru = size(Gu,1);
    Gbarpec = padarray(Gpec,[0 ((Nc-1)*nu+Nc*nx)], 0, 'post'); %zero-pad
    Gbar = [Gbar(1:nru,:); Gbarpec; Gbar(nru+1:end,:)];
    pbar = [pbar(1:nru,:); ppec; pbar(nru+1:end,:)];
end
%% Nonlinear solver
rel = version('-release');
rel = rel(1:4); %just the year
relnum = str2double(rel);
switch relnum
    case 2010
        options = optimset('Display', 'off');
    case 2011
        options = optimset(...
            'Algorithm', 'sqp', 'Display', 'off');
    case 2013
        options = optimoptions('fmincon', ...
            'Algorithm', 'sqp', 'Display', 'off');
    otherwise
        error('Can''t set solver options for this version of Matlab');
end
[Z ,FVAL,EXITFLAG, OUTPUT] = fmincon(@(z) 0.5*z'*Q_hat*z + q'*z, z0, ...
    Gbar, pbar, [], [], LB, UB, @nonlconfunc, options);
%% Return variables
X = reshape(Z, nu+nx,[]);
u = X(1:nu,:);
X = X(nu+1:end,:);
end