function aR = buildaR(omegas, sigmas, Gammas, thetas, R1N, C)
% Description: this function evaluates the residual linear accelerations
% needed for the computation of udot.

% Initialize aR vector
n = 5;      % n° of bodies
aR = zeros(3*n, 1);

% Retrieve Data from Input
R1N = reshape(R1N, 3, 3);
omega1 = omegas(1:3);
omega2 = omegas(4:6);
omega3 = omegas(7:9);
omega4 = omegas(10:12);
omega5 = omegas(13:15);
Gammas = reshape(Gammas, 3, 8);
GammaG1 = Gammas(:, 1);
GammaG2 = Gammas(:, 2);
GammaG3 = Gammas(:, 3);
GammaG4 = Gammas(:, 4);
GammadotG1 = Gammas(:, 5);
GammadotG2 = Gammas(:, 6);
GammadotG3 = Gammas(:, 7);
GammadotG4 = Gammas(:, 8);
theta1 = thetas(1);
theta2 = thetas(2);
theta3 = thetas(3);
theta4 = thetas(4);
sigma1 = sigmas(1);
sigma2 = sigmas(2);
sigma3 = sigmas(3);
sigma4 = sigmas(4);
R = C(1);
L = C(2);
lA = C(3);
lB = C(4);

% Define the vectors from each joint to CM of B1 in N f.o.r.
rG1_1 = R1N * [0 R 0]';
rG2_1 = R1N * [0 -R 0]';
rG3_1 = R1N * [0 0 -R]';
rG4_1 = R1N * [0 0 R]';

% Define the vectors from CM of B1 to each joint in N f.o.r.
r1_G1 = -rG1_1;
r1_G2 = -rG2_1;
r1_G3 = -rG3_1;
r1_G4 = -rG4_1;

% Define the vectors from CM of body j to the j-1th joint in N f.o.r.
r2_G1 = R1N * R21(theta1) * [0 -lA/2 0]';
r3_G2 = R1N * R31(theta2) * [0 -lA/2 0]';
r4_G3 = R1N * R41(theta3) * [0 0 -lB/2]';
r5_G4 = R1N * R51(theta4) * [0 0 -lB/2]';

% Build aR vector

aR(4:6) = skew(R1N*omega1)*skew(R1N*omega1)*r1_G1 ...
          - skew(R1N*omega1 + R1N*R21(theta1)*GammaG1*sigma1) * (skew(R1N*omega1 + R1N*R21(theta1)*GammaG1*sigma1) * r2_G1) ...
          - skew(R1N*R21(theta1) * (GammadotG1*sigma1 + skew(omega2)*GammaG1*sigma1))*r2_G1;

aR(7:9) = skew(R1N*omega1)*skew(R1N*omega1)*r1_G2 ...
          - skew(R1N*omega1 + R1N*R31(theta2)*GammaG2*sigma2) * (skew(R1N*omega1 + R1N*R31(theta2)*GammaG2*sigma2) * r3_G2) ...
          - skew(R1N*R31(theta2) * (GammadotG2*sigma2 + skew(omega3)*GammaG2*sigma2))*r3_G2;

aR(10:12) = skew(R1N*omega1)*skew(R1N*omega1)*r1_G3 ...
            - skew(R1N*omega1 + R1N*R41(theta3)*GammaG3*sigma3) * (skew(R1N*omega1 + R1N*R41(theta3)*GammaG3*sigma3) * r4_G3) ...
            - skew(R1N*R41(theta3) * (GammadotG3*sigma3 + skew(omega4)*GammaG3*sigma3))*r4_G3;

aR(13:15) = skew(R1N*omega1)*skew(R1N*omega1)*r1_G4 ...
            - skew(R1N*omega1 + R1N*R51(theta4)*GammaG4*sigma4) * (skew(R1N*omega1 + R1N*R51(theta4)*GammaG4*sigma4) * r5_G4) ...
            - skew(R1N*R51(theta4) * (GammadotG4*sigma4 + skew(omega5)*GammaG4*sigma4))*r5_G4;

 
end