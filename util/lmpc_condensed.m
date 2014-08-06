function [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_condensed(problem)
%LMPC_CONDENSED Compute the input sequence and predicted output using
%Linear MPC condensed (sequential) formulation.
%
%   [U, X, FVAL, EXITFLAG, OUTPUT] = LMPC_CONDENSED(PROBLEM). Compute the
%   input sequence using a sequential Linear MPC formulation over Nc
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
%   - uprev: the previously obtained input solution. This will be used to 
%   'warm start' the algorithm. This should be the u obtained at a previous 
%   step.
%   - Gx, px: constraint cell arrays or matrix for states. If Gx is a cell
%   array it must be of length Nc and comprised of nrxi-by-nx matrices and
%   px must be a Nc element cell array of nrxi element vectors. If Gx is a
%   matrix, is must be of size nrx-by-nx and px a nrx length vector. The
%   matrices of Gx describe constraints on combinations of states.
%       nrx - number of constraints on the combinations of states
%       nrxi - number of constraints on the combinations of states at
%       prediction i, where i is in [1, Nc].
%       Gx(i)*x(i) <= px(i)
%   - Gu, pu: constraint cell arrays or matrix for inputs. If Gu is a cell
%   array it must be of length Nc and comprised of nrui-by-nu matrices and
%   pu must be a Nc element cell array of nrui element vectors. If Gu is a
%   matrix, is must be of size nru-by-nu and pu a nru length vector. The
%   matrices of Gu describe constraints on combinations of states.
%       nru - number of constraints on the combinations of inputs
%       nrui - number of constraints on the combinations of inputs at
%       prediction i, where i is in [1, Nc].
%       Gu(i)*u(i) <= pu(i)
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
%% Argument processing
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
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
if isempty(uprev)
    %if there is no previous solution, use the reference as a start point
    uprev = uref;
else
    uprev = [uprev(:,1:end-1), uref(:,end)]; %shift the previous solution
end
%% QP definition
du(2,:) = -du(2,:); %convert negative constraints to positive
dx(2,:) = -dx(2,:); %convert negative constraints to positive
du = reshape(du',[],1); %reshape du into a vector
dx = reshape(dx',[],1); %reshape dx into a vector
Cx = [eye(nx); -eye(nx)];
Cu = [eye(nu); -eye(nu)];
%% ABbar
ABbarsmall = B;
ABbar = zeros(Nc*nx, Nc*nu);
for i = 1:Nc
    %Add another element to the block diagonal matrices
    %Add elements to ABbar
    for j = 1:(Nc + 1 - i)
        lstart = (i-1)*nx + (j-1)*nx + 1;
        lend = lstart + nx - 1;
        cstart = (j-1)*nu + 1;
        cend = cstart + nu - 1;
        ABbar(lstart:lend, cstart:cend) = ABbarsmall;
    end
    ABbarsmall = A*ABbarsmall;
end
%% Ap
Ap = [];
Apsmall = A;
for i = 1:Nc
    Ap = [Ap; Apsmall]; % Ap = [Ap A^i]
    Apsmall = Apsmall * A;
end
%% Gbar and pbar
Gubar = [];
pubar = [];
Gxbar = [];
pxbar = [];
if ~iscell(Gu)
    %transform to a Nc cell array
    Gu = repmat({Gu},[Nc 1]);
    pu = repmat({pu},[Nc 1]);
end
if ~iscell(Gx)
    %transform to a Nc cell array
    Gx = repmat({Gx},[Nc 1]); 
    px = repmat({px},[Nc 1]);
end
for i = 1:Nc
    %Append the next constraint, if it is empty there will be no change
    Gui = Gu{i};
    nrui = size(Gui, 1);
    if ~isempty(Gui)
        Guibar = [ zeros(nrui, (i-1)*nu), Gui, zeros(nrui, (Nc-i)*nu)];
    else
        Guibar = [];
    end
    Gubar = [Gubar; Guibar];
    pubar = [pubar; pu{i}];
    
    Gxi = Gx{i};
    nrxi = size(Gxi, 1);
    if ~isempty(Gxi)
        Gxibar = [ zeros(nrxi, (i-1)*nx), Gxi, zeros(nrxi, (Nc-i)*nx)];
    else
        Gxibar = [];
    end
    Gxbar = [Gxbar; Gxibar];
    pxbar = [pxbar; px{i}];
end
Gxbarhat = [];
if ~isempty(Gxbar)
    Gxbarhat = Gxbar * ABbar;
end
pxbarhat = [];
if ~isempty(pxbar)
    pxbarhat = pxbar - Gxbar*Ap*x0;
end
Gbar = [Gubar; Gxbarhat];
pbar = [pubar; pxbarhat];
%% Cz
Cxbar = [];
Cz2 = [];
for i = 1:Nc
    Cxbar = blkdiag(Cxbar, Cx);
    Cz2 = blkdiag(Cz2, Cu);
end
Cz1 = Cxbar*ABbar;
dxbar = repmat(dx,[Nc 1]);
dz1 = dxbar - Cxbar*Ap*x0;
dz2 = repmat(du,[Nc 1]);
C_hat = [Cz1; Cz2; Gbar];
d_hat = [dz1; dz2; pbar];
%% Qbar, Rbar
Qbar = [];
Rbar = [];
for i = 1:Nc
    Qbar = blkdiag(Qbar, Q);
    Rbar = blkdiag(Rbar, R);
end
%% Final
xrefbar = xref(:);
zrefbar = uref(:);
q = ABbar'*Qbar*Ap*x0 - ABbar'*Qbar*xrefbar - Rbar*zrefbar;
Q_hat = Rbar + ABbar'*Qbar*ABbar;
z0 = uprev(:);
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
            'Algorithm', 'active-set', 'Display', 'off');
    case 2013
        options = optimoptions('quadprog', ...
            'Algorithm', 'active-set', 'Display', 'off');
    otherwise
        error('Can''t set solver options for this version of Matlab');
end
[Z,FVAL,EXITFLAG, OUTPUT] = quadprog(Q_hat, q, C_hat, d_hat, [], [], ...
    [], [], z0, options);
%% Return variables
X = reshape(Z, nu,[]);
u = X(1:nu,:);
end
