function [C] = q2C(q0, q)
% Description: this function takes one attitude in input under the form of
% Quaternions and returns the rotation matrix using Directions Cosines.

C = (q0^2 - q'*q)*eye(3) + 2*(q*q') - 2*q0*skew(q);

end