clear
%% System
A = [ 1.0259     0.5040   0       0;
      1.0389     1.0259   0       0;
     -0.0006          0   1    0.05;
     -0.0247    -0.0006   0       1];
B = [-0.0013; -0.0504; 0.0006; 0.025];
x0 = [0.2; 0; 0; 0]; %Initial state
%% Cost matrices and constraints
Q = [ 1      0   0     0; 
      0   0.01   0     0; 
      0      0   1     0; 
      0      0   0  0.01];
R = 0.01;
dx = [0.2; inf; 1; inf;
      0.2; inf; 1; inf]; %state constraints, positive and -negative
du = [inf; inf]; %input constraints
%% Optimization problem definition
N = 100;
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
Cx = [eye(nx); -eye(nx)];
Cu = [eye(nu); -eye(nu)];
Csmall = blkdiag(Cu,Cx); 
Qsmall = blkdiag(R,Q);
Asmall = [-B eye(nx)];
bsmall= zeros(size(A,1),1);
C_hat = Csmall;
Q_hat = Qsmall;
A_hat = [-B eye(size(A,1))];
b_hat = A*x0;
for i = 1:N-1
    %Add another element to the block diagonal matrices
    C_hat = blkdiag(C_hat, Csmall);
    Q_hat = blkdiag(Q_hat, Qsmall);
    A_hat = blkdiag(A_hat, Asmall);
    %Add '-A' to the subdiagonal
    lines_l = i*nx + 1;
    lines_u = (i+1)*nx;
    cols_l = i*nu + (i-1)*nx + 1;
    cols_u = i*nu + i*nx;
    A_hat(lines_l: lines_u, cols_l: cols_u)= -A;
end
d_hat = repmat([du;dx], [N 1]);
b_hat = [b_hat; repmat(bsmall, [N-1 1])];
q = zeros(size(Q_hat,1),1);
options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', 'Display', 'final');
[Z,FVAL,EXITFLAG] = quadprog(Q_hat, q, C_hat, d_hat, A_hat, b_hat,[],[],[], options);
X = reshape(Z, 5,[]);
u = X(1,:);
X = X(2:5,:);
%% Plotting
t = 0:N-1;
figure(1);
subplot(4,2,1);
plot(t,u);
title('Input u');
grid on
subplot(4,2,2);
plot(t,X(1,:));
title('Arm position x_1');
grid on
subplot(4,2,4);
plot(t,X(2,:));
title('Arm speed x_2');
grid on
subplot(4,2,6);
plot(t,X(3,:));
title('Trolley position x_3');
grid on
subplot(4,2,8);
plot(t,X(4,:));
title('Trolley speed x_4');
grid on
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
titlestring = sprintf('Run with N=%d, Nc = %d',N,N);
text(0.5, 1, ['\bf ' titlestring],'HorizontalAlignment','center','VerticalAlignment', 'top')