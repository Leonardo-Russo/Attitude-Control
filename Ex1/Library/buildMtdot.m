function [Mtdot] = buildMtdot(xi, xidot)

% Description: this function creates the derivative in time of the Momenta
% Matrix Mt evaluating the derivative in time of each of its components.

DOF = 10;
md = 10;        % kg

JPdot = buildJPdot(md, xi, xidot);
SPdot = buildSPdot(md, xidot);

Mtdot = zeros(DOF, DOF);

Mtdot(4:6, 4:6) = JPdot;
Mtdot(4:6, 1:3) = skew(SPdot);
Mtdot(1:3, 4:6) = skew(-SPdot);

end