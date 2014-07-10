function [u, X, FVAL, EXITFLAG, OUTPUT] = lmpc_condensed(A, B, Q, R, Nc, du, dx, x0, xref)
%LMPC_CONDENSED Linear MPC in condensed form.
%   Explanation soon to come
%% Argument processing
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
if isempty(xref)
    xref = zeros(nx,1);
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
uref = zeros(nu,1); %should uref be 0?
xrefbar = repmat(xref, Nc, 1);
zrefbar = repmat(uref, Nc, 1);
q = ABbar'*Qbar*Ap*x0 - ABbar'*Qbar*xrefbar - Rbar*zrefbar;
Q_hat = Rbar + ABbar'*Qbar*ABbar;
%% QP solver
rel = version('-release');
rel = rel(1:4); %just the year
if strcmp(rel,'2011')
    options = optimset('Algorithm', 'active-set', 'Display', 'off'); % Matlab 2011
else
    options = optimoptions('quadprog', ...
        'Algorithm', 'active-set', 'Display', 'off'); %Matlab 2013
end
[Z,FVAL,EXITFLAG, OUTPUT] = quadprog(Q_hat, q, C_hat, d_hat, [], [], [], [], [], options);
%% Return variables
X = reshape(Z, nu,[]);
u = X(1:nu,:);
end

