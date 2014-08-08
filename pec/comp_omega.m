function [Omegatil_plus, Omegatil_minus] = comp_omega(U, N, P, rho0)
%COMP_OMEGA Compute Omega_tilde_plus matrix
%   Computes the Omega_tilde_plus matrix for time (k - 1).
m = size(U,1);
Omega = zeros(N*m, N*m);
Omegatil_minus = zeros(m*(N-1), m*(N-1));
for j = 1:P
    [Phi_minus, Phi] = comp_phi(U, N, j);
    Omega = Omega + Phi*Phi';
    Omegatil_minus = Omegatil_minus + Phi_minus*Phi_minus';
end
[~, Phi_kP] = comp_phi(U,N,P);
Omegatil_plus = Omega - Phi_kP * Phi_kP' - rho0 * eye(N*m);
Omegatil_minus = Omegatil_minus - rho0*eye(m*(N-1));
end
