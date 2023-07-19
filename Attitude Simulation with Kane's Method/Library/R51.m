function [R51] = R51(theta4)
% Description: this function returns the rotation matrix between the frames
% of reference from body 5 to body 1.

R51 = [-cos(theta4), 0, -sin(theta4);...
       0, 1, 0;...
       sin(theta4), 0, -cos(theta4)];

end