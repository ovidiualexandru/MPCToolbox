function x = quanser_disc_nl_euler(x0,u,h)
%QUANSER_DISC_NL_EULER Summary of this function goes here
%   Detailed explanation goes here
f = quanser_cont_nl([],[x0; u]);
xd = f(1:6);
x = x0 + h*xd;
end

