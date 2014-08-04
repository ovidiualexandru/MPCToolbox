function [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_condensed(A, B, Q, R, Nc,...
    du, dx, Gu, pu, Gx, px, x0, xref, uref, uprev)
%LMPC_CONDENSED Compute the input sequence and predicted output using
%Linear MPC condensed (sequential) formulation.
%
%   [U, X, FVAL, EXITFLAG, OUTPUT] = LMPC_CONDENSED(A, B, Q, R, Nc, DU, DX,
%   GU, PU, GX, PX, X0, XREF, UREF, UPREV). Compute the input sequence U 
%   for the model described by the model (A,B) using MPC condensed 
%   formulation, Nc as a control and prediction horizon, weighting matrices
%   Q and R. DU and DX describe the constraints for the inputs and states, 
%   X0 is the starting state, XREF is the state trajectory, UREF is the 
%   input trajectory and UPREV contains the last solution, used for a 'warm 
%   start'. The input sequence is returned in U, the predicted states in X.
%   FVAL, EXITFLAG and OUTPUT are returned by quadprog internally. See 
%   'help quadprog' for details.
%
%   Input arguments:
%   - A,B: the state-space matrices describing the system dynamic
%   - Q,R: the weighting matrices in the cost function
%   - Nc: the control horizon
%   - DU, DX: the constraint vectors for inputs and states. DU is a 2-by-nu
%   matrix containing constraints for inputs. First line is upper bound,
%   second is lower bound for each input. DX is a 2-by-nx matrix with
%   constraints for states. If the input/state has no lower bound, set it's
%   corresponding value to -Inf. Conversely, if the input/state has no
%   upper bound, set to Inf. 
%       nu - number of inputs, nx - number of states.
%   - GU, PU: constraint matrix and vector for inputs. GU is a nru-by-nu
%   matrix of combinations of constraints on the inputs and PU is a nru
%   vector. If there are no constraints, can be an empty matrix.
%       nru - number of constraints on the combinations of inputs
%       GU*u <= PU
%   - GX, PX: constraint matrix and vector for states. GX is a nrx-by-nx
%   matrix of combinations of constraints on the inputs and PU is a nrx
%   vector. If there are no constraints, can be an empty matrix.
%       nrx - number of constraints on the combinations of inputs
%       GX*x <= PX
%   - X0: the current( initial) state of the system
%   - XREF: the desired( reference) state. Must have nx lines, but can have
%   number of columns in the range [1, Nc].
%   - UREF: the reference input (stabilizing input). Must have nu lines,
%   but can have number of columns in the range [1, Nc]
%   - UPREV: the previously obtained input solution. This will be used as
%   a starting point for the algorithm. Can be an empty array if there 
%   is no previous solution, or the u obtained at a previous step.
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
Gxbar = [];
for i = 1:Nc
    Gubar = blkdiag(Gubar, Gu);
    Gxbar = blkdiag(Gxbar, Gx);
end
pubar = repmat(pu, [Nc 1]);
pxbar = repmat(px, [Nc 1]);
Gxbarhat = [];
if ~isempty(Gxbar)
    Gxbarhat = Gxbar * ABbar;
end
pxbarhat = [];
if ~isempty(pxbarhat)
    pxbarhat = pxbar - Ap*x0;
end
Gbar = blkdiag(Gubar, Gxbarhat);
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