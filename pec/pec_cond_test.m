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
det_val = zeros(Nv,1);
cond_val = zeros(Nv, 1);
for k = N+P:Nv
    Ux = Uk(:,1:k);
    Omegatil_plus = comp_omega(Ux, N, P, rho0);
    det_val(k) = det(Omegatil_plus);
    cond_val(k) = cond(Omegatil_plus);
end
%% Plotting
t = 1:Nv;
figure(1); plot(t,Uk(:,1:Nv)); title('Signal');grid on;
figure(2); plot(det_val); title('Determinant'); grid on;
figure(3); plot(cond_val); title('Condition number'); grid on;
