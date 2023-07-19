function [JP] = buildJP(xi)

% Description: this function builds the inertia Matrix JP from the position
% of the single damper mass xi.
% Note that all the constant quantities are defined inside this function
% and will not be present in the main program.

md = 10;        % kg

% Assume three equal Wheels
mw = 5;         % kg
Rw = 1;         % m
bw = 2;         % m

Is = mw * Rw^2 / 2;     % kg m^2
It = mw * Rw^2 / 4;     % kg m^2

% Define JB,P 
J1 = 350;       % kg m^2
J2 = 300;       % kg m^2
J3 = 400;       % kg m^2

JB_P = diag([J1, J2, J3]);

% Define JWi,wi for each wheel
JW3_w3_W3 = diag([It, It, Is]);
JW1_w1_W1 = diag([It, It, Is]);
JW2_w2_W2 = diag([It, Is, It]);

% Distance vectors between P and wi
rw1 = [0, 0, bw]';
rw2 = [0, bw, 0]';
rw3 = [0, 0, bw]';

% Distance vector between P and md
rd = [0, xi, 0]';

% Build Inertia Moments about point P
JW1_P_W1 = JW1_w1_W1 + mw *(dot(rw1, rw1)*eye(3) - rw1*rw1');
JW2_P_W2 = JW2_w2_W2 + mw *(dot(rw2, rw2)*eye(3) - rw2*rw2');
JW3_P_W3 = JW3_w3_W3 + mw *(dot(rw3, rw3)*eye(3) - rw3*rw3');

% Rotations to B from the f.o.r. of each wheel
alpha3 = pi/4;          % rad
alpha1 = 3*pi/4;        % rad

Rot_W3toB = [cos(alpha3), 0, sin(alpha3);...
              0, 1, 0;...
              -sin(alpha3), 0, cos(alpha3)];

Rot_W1toB = [cos(alpha1), 0, sin(alpha1);...
              0, 1, 0;...
              -sin(alpha1), 0, cos(alpha1)];

% Convert all Inertia Matrices into Body Frame
JW1_P = Rot_W1toB * JW1_P_W1 * Rot_W1toB';
JW2_P = JW2_P_W2;
JW3_P = Rot_W3toB * JW3_P_W3 * Rot_W3toB';

% Evaluate Inertia Matrix for the Damper
JD_P = md *(dot(rd, rd)*eye(3) - rd*rd');


% Finally sum all the terms
JP = JB_P + JW1_P + JW2_P + JW3_P + JD_P;

end