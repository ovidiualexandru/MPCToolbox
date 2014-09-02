%% Tunables
N = 10;
P = 15;
rho0 = .1;
%% Input signal
Uk = simout.U(:,450:600);
% Nv = 300;
Nv = size(Uk, 2);
% Uk = rand(1,Nv);
% Uk = Uk - mean(Uk); %remove dc component
% Uk = 0.2*Uk + 1.5; %scale
%% PEC check
% Ux = Uk(:,1:k);
[alpha, beta, gamma] = comp_pec(Uk, N, P, rho0);
