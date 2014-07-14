function [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_sparse(A, B, Q, R, Nc, du, dx, x0, xref)
%LMPC_SPARSE Calculate the input sequence and predicted output using Linear
%MPC sparse (simultaneous) formulation.
%   [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_sparse(A, B, Q, R, Nc, ...
%       du, dx, x0, xref). Calculate the inputs using MPC sparse
%       formulation.
%
%   Arguments:
%   - A,B: the state-space matrices describing the system dynamic
%   - Q,R: the weighting matrices in the cost function
%   - Nc: the control horizon
%   - du, dx: the constraint vectors for inputs and states. *du* is a 2-by-nu
%       matrix containing constraints for inputs. First line is lower
%       bound, second is upper bound for each input. *dx* is a 2-by-nx
%       matrix with constraints for states. If the input/state has no lower
%       bound, set it's corresponding value to -Inf. Conversely, if the
%       input/state has no upper bound, set to Inf. nu - number of inputs,
%       nx - number of states.
%   - x0: the current( initial) state of the system
%   - xref: the desired( reference) state
%   Output arguments:
%   - u: a nu-by-Nc matrix of computed inputs. u(:,1) must be used.
%   - X: a nx-by-Nc matrix of predicted states.
%   - FVAL: the object function value given by the numerical solver, 
%       quadprog.
%   - EXITFLAG: the exitflag from the solver. See 'help quadprog' for 
%       details. EXITFLAG is > 0 if a solution has been found.
%   - OUTPUT: the output from the solver. See 'help quadprog' for details.
%
%   Details for the sparse MPC formulation used can be found in 'Metode de
%       optimizare numerica'(romanian) by prof. I. Necoara, pg 237.

%% Argument processing
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
if isempty(xref)
    xref = zeros(nx,1);
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
uref = zeros(nu,1); %should uref be 0?
zsmall = [ uref; xref];
zref = repmat( zsmall, Nc,1);
q = -Q_hat*zref;
%% QP solver
rel = version('-release');
rel = rel(1:4); %just the year
if strcmp(rel,'2011')
    options = optimset('Algorithm', 'interior-point-convex', 'Display', 'off'); % Matlab 2011
else
    options = optimoptions('quadprog', ...
        'Algorithm', 'interior-point-convex', 'Display', 'off'); %Matlab 2013
end
[Z,FVAL,EXITFLAG, OUTPUT] = quadprog(Q_hat, q, [], [], A_hat, b_hat,LB,UB,[], options);
%% Return variables
X = reshape(Z, nu+nx,[]);
u = X(1:nu,:);
X = X(nu+1:end,:);