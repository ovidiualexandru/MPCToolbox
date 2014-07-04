function [u, X, FVAL, EXITFLAG, OUTPUT] = qp_fullstate(A, B, Q, R, Nc, du, dx, x0)
%% QP definition
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
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
q = zeros(size(Q_hat,1),1);
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