function dy = Dynamical_Model(t, y)

% Description: this function contains the Dynamical Model for the problem
% in question

global choice

n = length(y);
dy = zeros(n, 1);

% Define external forces and torques
Fext = zeros(3, 1);
Mp = zeros(3, 1);

nu = y(1:10);
xidot = nu(10);
xi = y(11);

JP = buildJP(xi);
SP = buildSP(xi);
Mt = buildMt(JP, SP);

B = buildB(nu, choice);
A = buildA(nu);
Mtdot = buildMtdot(xi, xidot);
c0 = buildc0(nu, Fext, Mp, xi, choice);

dy(1:10) = inv(Mt-B) * ((A*Mt - Mtdot)*nu + c0);
dy(11) = xidot;

end