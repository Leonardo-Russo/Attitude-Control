function [B] = buildB(nu, option)

% Description: this function builds the auxiliary vector c for the
% propagation.

% Import Known Quantities
md = 10;        % kg
Kw = 0.1;       % kg m^2/s

% Assume three equal Wheels
mw = 5;         % kg
Rw = 1;         % m
bw = 2;         % m

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

% Import Wheels Velocities from nu
omegas1 = nu(7);
omegas2 = nu(8);
omegas3 = nu(9);


% Create B Matrix
DOF = length(nu);
B = zeros(DOF, DOF);

switch option

    case '1'
        B(7, 4:6) = Is*a1';
        B(9, 4:6) = Is*a3';

    case '2'
        
end

end
