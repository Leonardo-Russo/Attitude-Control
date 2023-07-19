function [SP] = buildSP(xi)

% Description: this function builds the inertia Matrix JP from the position
% of the single damper mass xi.
% Note that all the constant quantities are defined inside this function
% and will not be present in the main program.

md = 10;        % kg

% Assume three equal Wheels
mw = 5;         % kg
Rw = 1;         % m
bw = 2;         % m

% Rotations to B from the f.o.r. of each wheel
alpha3 = pi/4;          % rad
alpha1 = 3*pi/4;        % rad

Rot_W3toB = [cos(alpha3), 0, sin(alpha3);...
              0, 1, 0;...
              -sin(alpha3), 0, cos(alpha3)];

Rot_W1toB = [cos(alpha1), 0, sin(alpha1);...
              0, 1, 0;...
              -sin(alpha1), 0, cos(alpha1)];


% Distance vectors between P and wi
rw1 = Rot_W1toB * [0, 0, bw]';
rw2 = [0, bw, 0]';
rw3 = Rot_W3toB * [0, 0, bw]';

% Distance vector between P and md
rd = [0, xi, 0]';


% Finally build Static Moment Vector
SP = mw*(rw1 + rw2 + rw3) + md*rd;


end