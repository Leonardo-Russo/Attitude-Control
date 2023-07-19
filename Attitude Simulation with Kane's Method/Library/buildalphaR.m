function alphaR = buildalphaR(omegas, sigmas, Gammas)
% Description: this function evaluates the angular residual accelerations starting from the
% analytical definition.

%%% I am considering alpha1R = 0

% Intialize alphaR vector
n = 5;      % n° of bodies
alphaR = zeros(3*n, 1);

% Retrieve Data from Input
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
sigma1 = sigmas(1);
sigma2 = sigmas(2);
sigma3 = sigmas(3);
sigma4 = sigmas(4);

% Build the alphaR vector

alphaR(4:6) = GammadotG1*sigma1 + skew(omega2)*GammaG1*sigma1;
alphaR(7:9) = GammadotG2*sigma2 + skew(omega3)*GammaG2*sigma2;
alphaR(10:12) = GammadotG3*sigma3 + skew(omega4)*GammaG3*sigma3;
alphaR(13:15) = GammadotG4*sigma4 + skew(omega5)*GammaG4*sigma4;


end