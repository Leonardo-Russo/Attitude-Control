function [V] = buildV(C, thetas, Gammas, R1N)
% Description: this functions build the 15x10 V matrix for the
% particular problem analysed.

% Initialize V Matrix
V = zeros(15, 10);

% Retrieve Data from Input
R = C(1);
L = C(2);
lA = C(3);
lB = C(4);
theta1 = thetas(1);
theta2 = thetas(2);
theta3 = thetas(3);
theta4 = thetas(4);
Gammas = reshape(Gammas, 3, 8);
GammaG1 = Gammas(:, 1);
GammaG2 = Gammas(:, 2);
GammaG3 = Gammas(:, 3);
GammaG4 = Gammas(:, 4);
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
r21 = r2_G1 + rG1_1;
r31 = r3_G2 + rG2_1;
r41 = r4_G3 + rG3_1;
r51 = r5_G4 + rG4_1;

% Build V Matrix

V(1:15, 8:10) = [eye(3); eye(3); eye(3); eye(3); eye(3)];

V(4:6, 1:3) = skew(r21) * R1N;
V(7:9, 1:3) = skew(r31) * R1N;
V(10:12, 1:3) = skew(r41) * R1N;
V(13:15, 1:3) = skew(r51) * R1N;

V(4:6, 4) = skew(r2_G1) * R1N * R21(theta1) * GammaG1;
V(7:9, 5) = skew(r3_G2) * R1N * R31(theta2) * GammaG2;
V(10:12, 6) = skew(r4_G3) * R1N * R41(theta3) * GammaG3;
V(13:15, 7) = skew(r5_G4) * R1N * R51(theta4) * GammaG4;

end