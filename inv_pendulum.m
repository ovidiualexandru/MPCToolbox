clear
%% System and restrictions
A = [ 1.0259    0.5040   0       0;
      1.0389    1.0259   0       0;
     -0.0006         0   1    0.05;
     -0.0247    -0.0006  0       1];
 
B = [-0.0013; -0.0504; 0.0006; 0.025];

x0 = [0.2; 0; 0; 0];
 
Q = [ 1      0   0     0; 
      0   0.01   0     0; 
      0      0   1     0; 
      0      0   0  0.01];
  
R = 1;

nu = size(B,2);
nx = size(A,1);

Cx =[ 1 0 0 0;
     -1 0 0 0;
    ];
Cu = [0; 0];
dx = [0.2; 0.2];
%% Construct C_hat and Q_hat
%Csmall = blkdiag(Cu,Cx);
Csmall = [zeros(2,1) Cx]; 
Qsmall = blkdiag(R,Q);
Asmall = [-B eye(size(A,1))];
bsmall= zeros(size(A,1),1);

N = 25;
C_hat = Csmall;
Q_hat = Qsmall;
A_hat = [-B eye(size(A,1))];
b_hat = A*x0;

for i = 2:N
    C_hat = blkdiag(C_hat, Csmall);
    Q_hat = blkdiag(Q_hat, Qsmall);
    A_hat = blkdiag(A_hat, Asmall);
    lines_l = nx + (i-2)*nx + 1;
    lines_u = lines_l + nx - 1;
    cols_l = nu + (i-2)*(nx + nu) + 1;
    cols_u = cols_l + nx - 1;
    %[lines_l lines_u; cols_l cols_u]
    A_hat(lines_l: lines_u, cols_l: cols_u)= -A;
end
dx_hat = repmat(dx, [N 1]);
b_hat = [b_hat; repmat(bsmall, [N-1 1])];
q = zeros(size(Q_hat,1),1);
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex');
[Z,FVAL,EXITFLAG] = quadprog(Q_hat, q, C_hat, dx_hat, A_hat, b_hat,[],[],[], options);
%% Plotting
figure(1);
plot(Z(1:5:end));
figure(2);
plot(Z(2:5:end));
figure(3);
plot(Z(3:5:end));
figure(4);
plot(Z(4:5:end));
figure(5);
plot(Z(5:5:end));