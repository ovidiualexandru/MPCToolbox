function [u, X, FVAL, EXITFLAG] = mpc_nonlin(handle_nlmodeld, h, Q, R, Nc, du, dx, x0)
%% Nonlinear constraints function
    function [C,Ceq] = nonlconfunc(z)
        %De for nu pot scapa
        % C e mereu 0
        % Ceq 'simuleaza' toate starile viitoare din X
        % ToDo:
        % - descompune X in [u x]
        C = zeros(size(z));
        Xl = reshape(z, nu+nx,[]);
        ul = Xl(1:nu,:);
        Xl = Xl(nu+1:end,:);
        for i = 1:Nc
            z = Xl(i,:);
            x = handle_nlmodeld(x, ul, h);
            Ceq(:) = z - x;
        end
    end
%% Argument processing
if ~isa(handle_nlmodeld, 'function_handle')
    error('handle_nlmodeld must be a function handle.');
end
%% QP definition
nu = size(du,2); %number of inputs
nx = size(dx,2); %number of states
du(2,:) = -du(2,:); %convert negative constraints to positive
dx(2,:) = -dx(2,:); %convert negative constraints to positive
du = reshape(du',[],1); %reshape du into a vector
dx = reshape(dx',[],1); %reshape dx into a vector
Cx = [eye(nx); -eye(nx)];
Cu = [eye(nu); -eye(nu)];
Csmall = blkdiag(Cu,Cx); 
Qsmall = blkdiag(R,Q);
C_hat = Csmall;
Q_hat = Qsmall;
for i = 1:Nc-1
    %Add another element to the block diagonal matrices
    C_hat = blkdiag(C_hat, Csmall);
    Q_hat = blkdiag(Q_hat, Qsmall);
end
d_hat = repmat([du;dx], [Nc 1]);
q = zeros(size(Q_hat,1),1);
%% Nonlinear solver
[Z ,FVAL,EXITFLAG] = fmincon(@(z) z'*Q_hat*z + q'*z, x0, C_hat, d_hat,[],[],[],[], @nonlconfunc);
%% Return variables
X = reshape(Z, nu+nx,[]);
u = X(1:nu,:);
X = X(nu+1:end,:);
end