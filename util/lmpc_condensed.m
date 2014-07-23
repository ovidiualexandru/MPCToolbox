function [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_condensed(A, B, Q, R, Nc,...
    du, dx, x0, xref, uref)
%LMPC_CONDENSED Compute the input sequence and predicted output using
%Linear MPC condensed (sequential) formulation.
%   [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_condensed(A, B, Q, R, Nc, ...
%       du, dx, x0, xref). Calculate the inputs using MPC condensed
%       formulation.
%
%   Arguments:
%   - A,B: the state-space matrices describing the system dynamic
%   - Q,R: the weighting matrices in the cost function
%   - Nc: the control horizon
%   - du, dx: the constraint vectors for inputs and states. *du* is a 
%       2-by-nu matrix containing constraints for inputs. First line is
%       upper bound, second is lower bound for each input. *dx* is a 
%       2-by-nx matrix with constraints for states. If the input/state has 
%       no lower bound, set it's corresponding value to -Inf. Conversely, 
%       if the input/state has no upper bound, set to Inf. 
%       nu - number of inputs, nx - number of states.
%   - x0: the current( initial) state of the system
%   - xref: the desired( reference) state. Must have nx lines, but can have
%       number of columns in the range [1, Nc].
%   - uref: the reference input (stabilizing input). Must have nu lines,
%       but can have number of columns in the range [1, Nc]
%   Output arguments:
%   - u: a nu-by-Nc matrix of computed inputs. u(:,1) must be used.
%   - X: a nu-by-Nc matrix of computed inputs, same as u.
%   - FVAL: the object function value given by the numerical solver, 
%       quadprog.
%   - EXITFLAG: the exitflag from the solver. See 'help quadprog' for 
%       details. EXITFLAG is > 0 if a solution has been found.
%   - OUTPUT: the output from the solver. See 'help quadprog' for details.
%
%   Details for the sparse MPC formulation used can be found in 'Metode de
%       optimizare numerica'(romanian) by prof. I. Necoara, pg 239.

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
C_hat = [Cz1; Cz2];
d_hat = [dz1; dz2];
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
    [], [], [], options);
%% Return variables
X = reshape(Z, nu,[]);
u = X(1:nu,:);
end