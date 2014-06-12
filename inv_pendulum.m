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
      -0.2; -inf; -inf; -inf]; %state constraints, positive and negative
du = [inf; -inf]; %input constraints
%% QP solve
N = 200;
Nc = 20;
nu = size(B,2); %number of inputs
nx = size(A,1); %number of states
X = zeros(nx, N); %save all states, for plotting
u = zeros(nu, N); %save all inputs
x = x0;
d = zeros(4,N);
dist_k = 60;
d(2, dist_k:end) = 0;
for i = 1:N
    X(:,i) = x; %save current state
    [ue, Xe] = qp_fullstate(A, B, Q, R, Nc, du, dx, x);    
    u(i) = ue(1); %use only the first command from predictions
    x = A*x + B*u(i); %compute next state
    x = x + 0.00.*rand(nx,1).*x + d(:,i); %add disturbance
end
%% Plotting
constraints_x = reshape(dx,[],2)';
constraints_u = du;

t = 0:N-1;
figure(1);
subplot(4,2,1);
plot(t,u);
rescaleYLim(gca, [constraints_u(2) constraints_u(1)]*1.1);
grid on
title('Input u');
line([0;N],[constraints_u(1);constraints_u(1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[constraints_u(2);constraints_u(2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,2);
plot(t,X(1,:));
rescaleYLim(gca, [constraints_x(2,1) constraints_x(1,1)]*1.1);
grid on
title('Arm position x_1');
line([0;N],[constraints_x(1,1);constraints_x(1,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[constraints_x(2,1);constraints_x(2,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,4);
plot(t,X(2,:));
rescaleYLim(gca, [constraints_x(2,2) constraints_x(1,2)]*1.1);
grid on
title('Arm speed x_2');
line([0;N],[constraints_x(1,2);constraints_x(1,2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[constraints_x(2,2);constraints_x(2,2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,6);
plot(t,X(3,:));
grid on
rescaleYLim(gca, [constraints_x(2,3) constraints_x(1,3)]*1.1); 
title('Trolley position x_3');
line([0;N],[constraints_x(1,3);constraints_x(1,3)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[constraints_x(2,3);constraints_x(2,3)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,8);
plot(t,X(4,:));
rescaleYLim(gca, [constraints_x(2,4) constraints_x(1,4)]*1.1);
grid on
title('Trolley speed x_4');
line([0;N],[constraints_x(1,4);constraints_x(1,4)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[constraints_x(2,4);constraints_x(2,4)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
line([dist_k;dist_k],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
titlestring = sprintf('Run with N=%d, Nc = %d',N,Nc);
text(0.5, 1, ['\bf ' titlestring],'HorizontalAlignment','center','VerticalAlignment', 'top')
