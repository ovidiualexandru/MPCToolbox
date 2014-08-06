function [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_sparse(problem)
%LMPC_SPARSE Compute the input sequence and predicted output using Linear
%MPC sparse (simultaneous) formulation.
%
%   [U, X, FVAL, EXITFLAG, OUTPUT] = LMPC_SPARSE(PROBLEM). Compute the
%   input sequence using a simultaneous Linear MPC formulation over Nc
%   samples for the linear system pair (A,B) with (Q,R) as weighting
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
%                             x(k+1) = A*x(k) + B*u(k)
%                             ul < u(k) < uu (input bounds)
%                             xl < x(k) < xu (state bounds)
%                             Gu*u(k) <= pu
%                             Gx*x(k) <= px, and
%                             du = [ul'; uu']
%                             dx = [xl'; xu']
%
%   Mandatory fields:
%   - A,B: the state-space matrices describing the system dynamic.
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
%
%   Output arguments:
%   - U: a nu-by-Nc matrix of computed inputs. u(:,1) must be used.
%   - X: a nu-by-Nc matrix of computed inputs, same as u.
%   - FVAL: the object function value given by the numerical solver,
%   quadprog.
%   - EXITFLAG: the exitflag from the solver. See 'help quadprog' for
%   details. EXITFLAG is > 0 if a solution has been found.
%   - OUTPUT: the output from the solver. See 'help quadprog' for details.
%
%   Details for the sparse MPC formulation used can be found in 'Metode de
%   optimizare numerica'(romanian) by prof. I. Necoara, pg 239.

%% Extract parameters from struct
% Required fields
A = problem.A;
B = problem.B;
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
%% Argument processing
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
if isempty(xref)
    xref = zeros(nx,1);
end
if isempty(uref)
    uref = zeros(nu,1);
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
Asmall = [-B eye(nx)];
bsmall= zeros(nx,1);
Q_hat = Qsmall;
A_hat = [-B eye(nx)];
b_hat = A*x0;
for i = 1:Nc-1
    %Add another element to the block diagonal matrices
    Q_hat = blkdiag(Q_hat, Qsmall);
    A_hat = blkdiag(A_hat, Asmall);
    %Add '-A' to the subdiagonal
    lines_l = i*nx + 1;
    lines_u = (i+1)*nx;
    cols_l = i*nu + (i-1)*nx + 1;
    cols_u = i*nu + i*nx;
    A_hat(lines_l: lines_u, cols_l: cols_u)= -A;
end
b_hat = [b_hat; repmat(bsmall, [Nc-1 1])];
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
%% QP solver
rel = version('-release');
rel = rel(1:4); %just the year
relnum = str2double(rel);
switch relnum
    case 2010
        options = optimset(...
            'Display', 'off', 'Diagnostics', 'off', 'LargeScale', 'off');
    case 2011
        options = optimset(...
            'Algorithm', 'interior-point-convex', 'Display', 'off');
    case 2013
        options = optimoptions('quadprog', ...
            'Algorithm', 'interior-point-convex', 'Display', 'off');
    otherwise
        error('Can''t set solver options for this version of Matlab');
end
[Z,FVAL,EXITFLAG, OUTPUT] = quadprog(Q_hat, q, Gbar, pbar, A_hat, b_hat, ...
    LB,UB, z0, options);
%% Return variables
X = reshape(Z, nu+nx,[]);
u = X(1:nu,:);
X = X(nu+1:end,:);