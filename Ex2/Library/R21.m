function [R21] = R21(theta1)
% Description: this function returns the rotation matrix between the frames
% of reference from body 2 to body 1.

R21 = [1, 0, 0;...
       0, -cos(theta1), sin(theta1);...
       0, -sin(theta1), -cos(theta1)];

end