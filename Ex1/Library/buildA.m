function [A] = buildA(nu)

% Description: this function builds the auxiliary matrix A for the
% propagation.

DOF = length(nu);
A = zeros(DOF, DOF);

vp = nu(1:3);
omega = nu(4:6);

A(1:3, 1:3) = -skew(omega);
A(4:6, 4:6) = -skew(omega);
A(4:6, 1:3) = -skew(vp);

end