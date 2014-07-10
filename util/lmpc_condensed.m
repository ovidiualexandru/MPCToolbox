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
Csmall = blkdiag(Cu,Cx);
for i = 1:Nc
    %Add another element to the block diagonal matrices
    %Add elements to ABbar
    ABbar = 0;
    Ap = 0; % Ap = [Ap A^i]
    Cxbar = 0; %Cxbar = blkdiag(Cxbar, Cx)
    Cz2 = 0; %Cz2 = blkdiag(Cz2, Cu)
end
Cz1 = 0; %Cxbar*ABbar
dz1 = 0; %dxbar - Cxbar*Ap*x0
dz2 = 0; %repmat(du,[Nc 1]);
C_hat = 0; % [Cz1; Cz2]
d_hat = 0; % [dz1; dz2]
uref = zeros(nu,1); %should uref be 0?
xrefbar = repmat(xref, Nc, 1);
zrefbar = repmat(uref, Nc, 1);
q = 0; % ABbar*Qbar*Ap*x0 - ABbar'*Qbar*xrefbar - Rbar*zrefbar
Qhat = 0; %Rbar + ABbar'*Qbar*ABbar;
%% QP solver
rel = version('-release');
rel = rel(1:4); %just the year
if strcmp(rel,'2011')
    options = optimset('Algorithm', 'interior-point-convex', 'Display', 'off'); % Matlab 2011
else
    options = optimoptions('quadprog', ...
        'Algorithm', 'interior-point-convex', 'Display', 'off'); %Matlab 2013
end
[Z,FVAL,EXITFLAG, OUTPUT] = quadprog(Q_hat, q, C_hat, d_hat, [], [], [], [], [], options);
%% Return variables
X = reshape(Z, nu+nx,[]);
u = X(1:nu,:);
X = X(nu+1:end,:);
end

