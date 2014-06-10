clear
%% System and restrictions
A = [ 1.0259     0.5040   0       0;
      1.0389     1.0259   0       0;
     -0.0006          0   1    0.05;
     -0.0247    -0.0006   0       1];
 
B = [-0.0013; -0.0504; 0.0006; 0.025];

x0 = [0.2; 0; 0; 0];
 
Q = [ 1      0   0     0; 
      0   0.01   0     0; 
      0      0   1     0; 
      0      0   0  0.01];
  
R = 10;

nu = size(B,2); %number of inputs
nx = size(A,1); %number of states

Cx =[ 1 0 0 0;
     -1 0 0 0];
Cu = [0; 0];
dx = [0.2; 0.2];
%% Construct C_hat and Q_hat
Csmall = [zeros(2,1) Cx]; 
Qsmall = blkdiag(R,Q);
Asmall = [-B eye(size(A,1))];
bsmall= zeros(size(A,1),1);

N = 25;
C_hat = Csmall;
Q_hat = Qsmall;
A_hat = [-B eye(size(A,1))];
b_hat = A*x0;

for i = 1:N-1
    C_hat = blkdiag(C_hat, Csmall);
    Q_hat = blkdiag(Q_hat, Qsmall);
    A_hat = blkdiag(A_hat, Asmall);
    lines_l = i*nx + 1;
    lines_u = (i+1)*nx;
    cols_l = i*nu + (i-1)*nx + 1;
    cols_u = i*nu + i*nx;
    %[lines_l lines_u; cols_l cols_u]
    A_hat(lines_l: lines_u, cols_l: cols_u)= -A;
end
dx_hat = repmat(dx, [N 1]);
b_hat = [b_hat; repmat(bsmall, [N-1 1])];
q = zeros(size(Q_hat,1),1);
options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', 'Display', 'off');
[Z,FVAL,EXITFLAG] = quadprog(Q_hat, q, C_hat, dx_hat, A_hat, b_hat,[],[],[], options);
%% Plotting
X = reshape(Z, 5,[]);
X = X';
u = X(:,1);
X = X(:,2:5);
figure(1);
subplot(4,2,1);
plot(u);
title('Input u');
grid on
subplot(4,2,2);
plot(X(:,1));
title('Arm position x_1');
grid on
subplot(4,2,4);
plot(X(:,2));
title('Arm speed x_2');
grid on
subplot(4,2,6);
plot(X(:,3));
title('Trolley position x_3');
grid on
subplot(4,2,8);
plot(X(:,4));
title('Trolley speed x_4');
grid on