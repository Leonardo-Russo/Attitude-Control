function dX = DynamicalModel(t, X)
% Description: this function contains the Dynamical Model used for the
% propagation in time of this problem.

% Initialize Derivative Value
dX = zeros(length(X), 1);

% Recall Known Quantities
global mu J Is sign_qe0_0 section omega_dot_max

% Import Data from Input
w = X(1:3);
omegas = X(4:7);
r_vect = X(8:10);
v_vect = X(11:13);
xe = X(14:17);

r = norm(r_vect);
qe0 = xe(1);
qe = xe(2:4);

% Evaluate Angular Momentum h
h_vect = cross(r_vect, v_vect);
h = norm(h_vect);

% Gain Parameters
omega_n = 0.03;     % rad/s
xi = 1;
c1 = 2 * omega_n^2;
c2 = 2*xi*omega_n/c1;
invA = c1 * eye(3);
B = c2 * eye(3);

% Define Rotation Matrix from ECI to Perifocal f.o.r.
COE = rvECI2coe(r_vect, v_vect, mu);
RNP = RotECI2Peri(COE);

% Evaluate Radial Velocity
v_vectP = RNP * v_vect;
vr = v_vectP(1);

% Define the Commanded Omega and OmegaDot
wc = [0 -h/r^2 0]';
wc_dot = [0 2*h/r^3*vr 0]';

% Define Desired Omega
RCB = q2C(qe0, qe);
we = w - RCB * wc;
wd = w - wc;


% Define aj as the axis of each RW
a1 = [1 0 0]';
a2 = [0 1 0]';
a3 = [0 0 1]';
a4 = [sqrt(3)/3 sqrt(3)/3 sqrt(3)/3]';
as = [a1, a2, a3, a4];

% Build Necessary Matrices
Ja = buildJa(Is, as);

% Define Mc
Mc = 0;

% Compute Commanded Torque
Tc = skew(w)*J*w - Mc + J*wc_dot - J*invA*B*wd - sign_qe0_0*J*invA*qe;

% Assign Derivative Values
dX(8:10) = v_vect;
dX(11:13) = -mu/r^3 * r_vect;
dX(14) = -0.5 * qe' * wd;
dX(15:17) = 0.5 * (qe0*eye(3) + skew(qe)) * we;


if section == 'A' | section == 'B'
    dX(4:7) = -Ja' * inv(Ja * Ja') * Tc;
elseif section == 'C'
    upd = -Ja' * inv(Ja * Ja') * Tc;
    for j = 4 : 7
        if upd(j-3) > omega_dot_max || upd(j-3) < -omega_dot_max
            dX(j) = omega_dot_max * sign(upd(j-3));
        else
            dX(j) = upd(j-3);
        end
    end
end

dX(1:3) = inv(J) * (Mc - skew(w)*J*w - skew(w)*Ja*omegas - Ja*dX(4:7));


end