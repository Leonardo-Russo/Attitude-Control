function F = buildF(C, thetas, rp1, R1N, M, m)
% Description: this function evaluates the Forces applied on the bodies

% Import Known Quantities
mu = 398600.4415 * 1e9;     % m^3/s^2

% Retrieve Data from Input
R = C(1);
L = C(2);
lA = C(3);
lB = C(4);
theta1 = thetas(1);
theta2 = thetas(2);
theta3 = thetas(3);
theta4 = thetas(4);
R1N = reshape(R1N, 3, 3);

% Define the vectors from each joint to CM of B1 in N f.o.r.
rG1_1 = R1N * [0 R 0]';
rG2_1 = R1N * [0 -R 0]';
rG3_1 = R1N * [0 0 -R]';
rG4_1 = R1N * [0 0 R]';

% Define the vectors from CM of body j to the j-1th joint in N f.o.r.
r2_G1 = R1N * R21(theta1) * [0 -lA/2 0]';
r3_G2 = R1N * R31(theta2) * [0 -lA/2 0]';
r4_G3 = R1N * R41(theta3) * [0 0 -lB/2]';
r5_G4 = R1N * R51(theta4) * [0 0 -lB/2]';

% Define path vectors in N f.o.r.
r12 = -(r2_G1 + rG1_1);
r13 = -(r3_G2 + rG2_1);
r14 = -(r4_G3 + rG3_1);
r15 = -(r5_G4 + rG4_1);

% Define Inertial Path Vectors
r2_vect = rp1 + r12;
r3_vect = rp1 + r13;
r4_vect = rp1 + r14;
r5_vect = rp1 + r15;

r1 = norm(rp1);
r2 = norm(r2_vect);
r3 = norm(r3_vect);
r4 = norm(r4_vect);
r5 = norm(r5_vect);

% Compute Gravitational Forces
F1 = -mu*M/r1^3 * rp1;
F2 = -mu*m/r2^3 * r2_vect;
F3 = -mu*m/r3^3 * r3_vect;
F4 = -mu*m/r4^3 * r4_vect;
F5 = -mu*m/r5^3 * r5_vect;

F = [F1; F2; F3; F4; F5];

end