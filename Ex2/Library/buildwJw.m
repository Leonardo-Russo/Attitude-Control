function wJw = buildwJw(omegas, J)
% Description: this function creates the support vector wJw needed in the
% computation of udot.

% Initialize wJw vector
n = 5;      % n° of bodies
wJw = zeros(3*n, 1);

% Retrieve Data from Input
omega1 = omegas(1:3);
omega2 = omegas(4:6);
omega3 = omegas(7:9);
omega4 = omegas(10:12);
omega5 = omegas(13:15);
J = reshape(J, 15, 15);
J1 = J(1:3, 1:3);
J2 = J(4:6, 4:6);
J3 = J(7:9, 7:9);
J4 = J(10:12, 10:12);
J5 = J(13:15, 13:15);

% Build the wJw vector

wJw(1:3) = skew(omega1)*J1*omega1;
wJw(4:6) = skew(omega2)*J2*omega2;
wJw(7:9) = skew(omega3)*J3*omega3;
wJw(10:12) = skew(omega4)*J4*omega4;
wJw(13:15) = skew(omega5)*J5*omega5;


end