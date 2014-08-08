N = 10;
P = 15;
rho0 = 0.1;

Nv = size(U,2);
Nv = 300;
Uk = 10*(rand(1,Nv) - 0.5);
det_val = zeros(Nv,1);
cond_val = zeros(Nv, 1);
for k = N+P:Nv
    Ux = Uk(:,1:k);
    Omegatil_plus = comp_omega(Ux, N, P, rho0);
    det_val(k) = det(Omegatil_plus);
    cond_val(k) = cond(Omegatil_plus);
end
t = 1:Nv;
figure(1); plot(t,Uk(1:Nv)); grid on;
figure(2); plot(det_val); grid on;
figure(3); plot(cond_val); grid on;
