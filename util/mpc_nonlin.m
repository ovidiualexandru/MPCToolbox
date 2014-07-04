function [u, X, FVAL, EXITFLAG, OUTPUT] = mpc_nonlin(handle_nlmodeld, h, Q, R, Nc, du, dx, x0)
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
            x = handle_nlmodeld(x, ul(:,i), h);
            dif = b - x; %x_{k+1} - f(x_k)
            ceq = [zeros(nu, 1); dif]; %append zeros because inputs have no eq constraints
            idx_start = (i-1)*(nx+nu)+1;
            idx_end = i*(nx+nu);
            Ceq(idx_start:idx_end) = ceq;
        end
    end
%% Argument processing
if ~isa(handle_nlmodeld, 'function_handle')
    error('handle_nlmodeld must be a function handle.');
end
%% QP definition
nu = size(du,2); %number of inputs
nx = size(dx,2); %number of states
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
q = zeros(size(Q_hat,1),1);
z0 = zeros(size(Q_hat,1),1);
%% Nonlinear solver
options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', 'Display', 'off'); %Matlab 2013
[Z ,FVAL,EXITFLAG, OUTPUT] = fmincon(@(z) z'*Q_hat*z + q'*z, z0, [], [],[],[],LB,UB, @nonlconfunc, options);
%% Return variables
X = reshape(Z, nu+nx,[]);
u = X(1:nu,:);
X = X(nu+1:end,:);
end