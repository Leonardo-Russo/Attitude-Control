function [lambda, phi] = rect2polar(r_vect)
% Description:
% Transformation from cartesian coordinates to polar.

% Note:
% All angles are transformed into [deg].

% Obtain values from input
x = r_vect(1);
y = r_vect(2);
z = r_vect(3);

% Compute output values
r = norm(r_vect);
phi = (asin(z/r)) * (180/pi());         % deg
lambda = atan2(y, x) * (180/pi());      % deg
