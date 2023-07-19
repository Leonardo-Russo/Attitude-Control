function [q0, q] = C2q(C)
% Description: this function obtains the quaternions from the Direction
% Cosines matrix provided as input.

q0 = 0.5 * sqrt(1 + C(1,1) + C(2,2) + C(3,3));

q = 1/(4*q0) * [C(2,3)-C(3,2); C(3,1)-C(1,3); C(1,2)-C(2,1)];

end