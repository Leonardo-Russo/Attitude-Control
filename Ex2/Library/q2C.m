function [C] = q2C(q0, q)
% Description: this function takes one attitude in input under the form of
% Quaternions and returns the rotation matrix using Directions Cosines.

q1 = q(1);
q2 = q(2);
q3 = q(3);

C = zeros(3, 3);

C(1, 1) = 1 - 2*(q2^2+q3^2);
C(1, 2) = 2*(q1*q2+q0*q3);
C(1, 3) = 2*(q1*q3-q0*q2);
C(2, 1) = 2*(q2*q1-q0*q3);
C(2, 2) = 1 - 2*(q1^2+q3^2);
C(2, 3) = 2*(q2*q3+q0*q1);
C(3, 1) = 2*(q3*q1+q0*q2);
C(3, 2) = 2*(q3*q2-q0*q1);
C(3, 3) = 1 - 2*(q1^2+q2^2);

end