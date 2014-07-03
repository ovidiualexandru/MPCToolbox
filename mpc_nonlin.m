% Objective function for fmincon

function [u,X, FVAL, EXITFLAG] = mpc_nonlin(x0)
%Cum specific x0?
    X = fmincon(@objfunc, x0, A,B,[],[],[],[], @nonlconfunc);
end

function F = objfunc(X)
% X has concatenations L of [u;x]
%de for pot scapa daca fac Q_hat = blkdiag(R,Q)
F = 0;
for i = 1:L
    F = F + x'*Q*x + u'*R*u;
end
end

%Nonlinear constraints for fmincon
function [C,Ceq] = nonlconfunc(X)
%De for nu pot scapa
% C e mereu 0
% Ceq 'simuleaza' toate starile viitoare din X
% ToDo:
% - descompune X in [u x]
C = 0;
for i = 1:L
    x = nonlindisc(x,u,h);
    Ceq(:) = X - x;
end
end

function xr = nonlindisc(x,u,h)
[Tout, Yout] = ode45(@quanser_cont_nl, [0 h], [x; u]); %f(xk, uk)
xr = Yout(end, 1:6)'; %get new state, i.e. x = x(k)
end


