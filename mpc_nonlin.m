% Objective function for fmincon

function [u,Z, FVAL, EXITFLAG] = mpc_nonlin(x0, handle_nonlindisc)
%Nonlinear constraints for fmincon
    function [C,Ceq] = nonlconfunc(z)
        %De for nu pot scapa
        % C e mereu 0
        % Ceq 'simuleaza' toate starile viitoare din X
        % ToDo:
        % - descompune X in [u x]
        C = 0;
        for i = 1:L
            x = nonlindisc(x,u,h);
            Ceq(:) = z - x;
        end
    end

if ~isa(handle_nonlindisc, 'function_handle')
    error('No function handle passed');
end
%% Solve nonlinear problem
Z = fmincon(@(z) z'*Q_hat*z, x0, A,B,[],[],[],[], @nonlconfunc);
end




