function [R31] = R31(theta2)
% Description: this function returns the rotation matrix between the frames
% of reference from body 3 to body 1.

R31 = [1, 0, 0;...
       0, cos(theta2), -sin(theta2);...
       0, sin(theta2), cos(theta2)];

end