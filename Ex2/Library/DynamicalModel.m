function dy = DynamicalModel(t, y)
% Description: this function contains the evolution in time of all the
% state variables of the problem.

% Recall Global variables
global mB1 mBj R L r lA lB K1 K2

% From the analysis of motion of each joint and the connected body
GammaG1 = [1 0 0]';
GammaG2 = [1 0 0]';
GammaG3 = [0 1 0]';
GammaG4 = [0 1 0]';
GammaG1dot = zeros(3, 1);
GammaG2dot = zeros(3, 1);
GammaG3dot = zeros(3, 1);
GammaG4dot = zeros(3, 1);
Gammas = [GammaG1, GammaG2, GammaG3, GammaG4, GammaG1dot, GammaG2dot, GammaG3dot, GammaG4dot];

% Import Data from Input State
omega1 = y(1:3);
sigmas = y(4:7);
vp1 = y(8:10);
rp1 = y(11:13);
thetas = y(14:17);
q0 = y(18);
q = y(19:21);
u = y(1:10);

% Compute Rotation matrix from 1 to N
R1N = q2C(q0, q)';
R1N = reshape(R1N, 9, 1);

% Initialize the Derivative values
dy = zeros(length(y), 1);

% Compute Necessary Quantities
M = buildM(mB1, mBj);
J = buildJ(mB1, mBj, R, L, r, lA, lB);
C = [R, L, lA, lB]';

T = buildT(sigmas, thetas, Gammas, K1, K2);
F = buildF(C, thetas, rp1, R1N, mB1, mBj);

Omega = buildOmega(thetas, Gammas);
V = buildV(C, thetas, Gammas, R1N);

omegas = Omega * u;

alphaR = buildalphaR(omegas, sigmas, Gammas);
aR = buildaR(omegas, sigmas, Gammas, thetas, R1N, C);

Jvec = reshape(J, 15*15, 1);    % reshaped vector needed for wJw(...)
wJw = buildwJw(omegas, Jvec);

udot = (Omega' * J * Omega + V' * M * V) \ (Omega'*(T - J*alphaR - wJw) + V'*(F - M*aR));

% Fill the Derivative Vector
dy(1:10) = udot;
dy(11:13) = vp1;
dy(14:17) = sigmas;
dy(18) = -0.5 * dot(q, omega1);
dy(19:21) = 0.5 * (q0*omega1 + cross(q, omega1));

end