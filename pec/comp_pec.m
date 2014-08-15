function [alpha, beta, gamma] = comp_pec(U, N, P, rho0)
%COMP_PEC Compute the PEC condition terms
%   [alfa, beta, gamma] = COMP_PEC(U, N, P, rho0)
m = size(U,1);
%% Alpha
[~, Omegatil_minus] = comp_omega(U, N, P, rho0);
Phi_minus = comp_phi(U, N, 1);
invOmegatil_minus = inv(Omegatil_minus);
alpha = 1 - Phi_minus'*invOmegatil_minus*Phi_minus;
%% Beta
beta = 0;
for j = 1:P-1
    uj = U(:,end-(j-1));
    Phi_minusj = comp_phi(U,N,j);
    beta = beta - uj*Phi_minusj'*invOmegatil_minus*Phi_minus;
end
%% Gamma
sum1 = 0;
for j = 1:P-1
    uj = U(:,end-(j-1));
    sum1 = sum1 + uj*uj';
end
sum1 = sum1 - rho0*eye(m);

sum2 = 0;
for j = 1:P-1
    uj = U(:,end-(j-1));
    Phi_minusj = comp_phi(U,N,j);
    sum2 = sum2 + uj*Phi_minusj';
end

sum3 = 0;
for j= 1:P-1
    uj = U(:,end-(j-1));
    Phi_minusj = comp_phi(U,N,j);
    sum3 = sum3 + Phi_minusj*uj';
end
gamma = sum1 - sum2 * invOmegatil_minus * sum3;
end
