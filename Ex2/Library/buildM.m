function [M] = buildM(mB1, mBj)
% Description: this function generates the M matrix from the individual
% masses of the bodies.

n = 5;      % n° of bodies

M = zeros(3*n, 3*n);

M(1:3, 1:3) = mB1 * eye(3);
M(4:6, 4:6) = mBj * eye(3);
M(7:9, 7:9) = mBj * eye(3);
M(10:12, 10:12) = mBj * eye(3);
M(13:15, 13:15) = mBj * eye(3);

end