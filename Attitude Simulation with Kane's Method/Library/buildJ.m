function [J] = buildJ(mB1, mBj, R, L, r, lA, lB)
% Description: this function creates the J matrix from the geometrical
% parameters and mass of each body.

% Initialize J Matrix
n = 5;      % n° of bodies
J = zeros(3*n, 3*n);

% Compute Inertia Moments
Ia1 = mB1 * R^2 / 2;
It1 = 1/12 * mB1 * (L^2 + 3*R^2);
Iaj = mBj * r^2 / 2;
ItjA = 1/12 * mBj * (lA^2 + 3*r^2);
ItjB = 1/12 * mBj * (lB^2 + 3*r^2);

% Build J Matrix
J(1:3, 1:3) = diag([Ia1, It1, It1]);
J(4:6, 4:6) = diag([ItjA, Iaj, ItjA]);
J(7:9, 7:9) = diag([ItjA, Iaj, ItjA]);
J(10:12, 10:12) = diag([ItjB, ItjB, Iaj]);
J(13:15, 13:15) = diag([ItjB, ItjB, Iaj]);

end