function [Mt] = buildMt(JP, SP)

% Descriptio: this function creates the Momenta Matrix Mt with the values
% provided at each integration step.
% Input:

DOF = 10;
Mt = zeros(DOF, DOF);

mb = 100;       % kg
md = 10;        % kg

% Assume three equal Wheels
mw = 5;         % kg
Rw = 1;         % m
bw = 2;         % m

M = mb + 3*mw + md;     % total mass M

Is = mw * Rw^2 / 2;     % kg m^2
It = mw * Rw^2 / 4;     % kg m^2

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

% Define the versors for wheels and dampers
a1 = rw1/norm(rw1);
a2 = rw2/norm(rw2);
a3 = rw3/norm(rw3);
n = [0 1 0]';


% Finally, build the Mt matrix
Mt(1:3, 1:3) = M * eye(3);
Mt(4:6, 4:6) = JP;
Mt(4:6, 1:3) = skew(SP);
Mt(1:3, 4:6) = skew(-SP);
Mt(7, 4:6) = Is*a1';
Mt(8, 4:6) = Is*a2';
Mt(9, 4:6) = Is*a3';
Mt(4:6, 7) = Is*a1;
Mt(4:6, 8) = Is*a2;
Mt(4:6, 9) = Is*a3;
Mt(7:10, 7:10) = diag([Is Is Is md]);
Mt(10, 1:3) = md*n';
Mt(10, 4:6) = 0;        % because skew(b) is null
Mt(1:3, 10) = md*n;
Mt(4:6, 10) = 0;        % because skew(b) is null


end