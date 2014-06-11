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
dx = [0.2; inf; inf; inf;
      0.2; inf; inf; inf]; %state constraints, positive and -negative
du = [inf; inf]; %input constraints
%% QP solve
N = 200;
Nc = 20;
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
X = zeros(nx, N); %save all states, for plotting
u = zeros(nu, N); %save all inputs
x = x0;
for i = 1:N
    %fprintf('Iteration: %d\n',i);
    X(:,i) = x; %save current state
    [ue, Xe] = qp_fullstate(A, B, Q, R, Nc, du, dx, x);    
    u(i) = ue(1); %use only the first command from predictions
    x = A*x + B*u(i); %compute next state
    x = x + 0.01.*rand(nx,1).*x; %add disturbance
end
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