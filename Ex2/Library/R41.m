function [R41] = R41(theta3)
% Description: this function returns the rotation matrix between the frames
% of reference from body 4 to body 1.

R41 = [cos(theta3), 0, sin(theta3);...
       0, 1, 0;...
       -sin(theta3), 0, cos(theta3)];

end