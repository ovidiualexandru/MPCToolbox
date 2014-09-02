function [Phi_minus, Phi] = comp_phi(U, N, j)
%PHI_MINUS Phi vector for time instance (k - j)
%   [Phi_minus, Phi] = COMP_PHI(U, N, j) returns the vectors Phi_minus and
%   Phi for time instance (k - j), using input values from U and N as 
%   number of variables that should be estimated.
%
%   Input arguments:
%   - U: Input value matrix. Is a m-by-M matrix with input values ordered
%   by column. <m> is the number of inputs. M is a number greater than
%   (j + N - 1). It is assumed that the last column in U correspons to time
%   (k - 1).
%   - N: excitation order
%   - j: number of samples to go back. The sample matrix is assumed to end
%   at time (k-1) to j must be >= 1.
%
%   Output arguments:
%   - Phi_minus: the phi_minus vector corresponding to time (k-j).
%   - Phi_minus: the phi vector corresponding to time (k-j).
%
m = size(U,1);
% Remember: end = k-1
U = U(:, (end-j-N+2):(end-j+1));
Phi = (U(:,end:-1:1));
Phi = Phi(:);
Phi_minus = Phi(1:end-m);
end
