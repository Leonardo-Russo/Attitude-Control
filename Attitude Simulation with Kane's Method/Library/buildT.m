function T = buildT(sigmas, thetas, Gammas, K1, K2)
% Description: this function evaluates the Torques on each joint and
% returns them in a 15x1 vector.

% Retrieve Data from Input
theta1 = thetas(1);
theta2 = thetas(2);
theta3 = thetas(3);
theta4 = thetas(4);
sigma1 = sigmas(1);
sigma2 = sigmas(2);
sigma3 = sigmas(3);
sigma4 = sigmas(4);
Gammas = reshape(Gammas, 3, 8);
GammaG1 = Gammas(:, 1);
GammaG2 = Gammas(:, 2);
GammaG3 = Gammas(:, 3);
GammaG4 = Gammas(:, 4);

% Compute each Torque
T2 = (-K1*theta1 - K2*sigma1) * GammaG1;
T3 = (-K1*theta2 - K2*sigma2) * GammaG2;
T4 = (-K1*theta3 - K2*sigma3) * GammaG3;
T5 = (-K1*theta4 - K2*sigma4) * GammaG4;
T1 = -(R21(theta1)*T2 + R31(theta2)*T3 + R41(theta3)*T4 + R51(theta4)*T5);

T = [T1; T2; T3; T4; T5];

end