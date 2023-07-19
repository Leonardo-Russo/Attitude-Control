function [Omega] = buildOmega(thetas, Gammas)
% Description: this functions build the 15x10 Omega matrix for the
% particular problem analysed.

% Initialize Omega Matrix
Omega = zeros(15, 10);

% Retrieve Data from input
theta1 = thetas(1);
theta2 = thetas(2);
theta3 = thetas(3);
theta4 = thetas(4);
Gammas = reshape(Gammas, 3, 8);
GammaG1 = Gammas(:, 1);
GammaG2 = Gammas(:, 2);
GammaG3 = Gammas(:, 3);
GammaG4 = Gammas(:, 4);

% Build Omega Matrix

Omega(1:3, 1:3) = eye(3);

Omega(4:6, 1:3) = R21(theta1)';
Omega(7:9, 1:3) = R31(theta2)';
Omega(10:12, 1:3) = R41(theta3)';
Omega(13:15, 1:3) = R51(theta4)';

Omega(4:6, 4) = GammaG1;
Omega(7:9, 5) = GammaG2;
Omega(10:12, 6) = GammaG3;
Omega(13:15, 7) = GammaG4;

end