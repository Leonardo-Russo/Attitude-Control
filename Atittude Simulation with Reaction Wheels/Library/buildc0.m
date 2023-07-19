function c0 = buildc0(nu, Fext, Mp, xi, option)

% Description: this function builds the auxiliary vector c for the
% propagation.

% Recall known Quantities
md = 10;        % kg
n = [0 1 0]';   % orientation of the damper
Kd = 5.5;       % kg/s^2
Cd = 30;        % kg/s
Kw = 0.1;       % kg m^2/s

DOF = length(nu);
c0 = zeros(DOF, 1);

vp = nu(1:3);
omega = nu(4:6);
omegas = nu(7:9);
xidot = nu(10);

rd = xi * n;    % position of the damper mass

c0(1:3) = Fext;
c0(4:6) = Mp;
c0(10) = md * omega' * skew(n) * (vp - skew(rd)*omega) - Cd*xidot -Kd*xi;

switch option
    
    case '1'
        c0(8) = -Kw*omegas(2);
    
    case '2'
        c0(7:9) = -Kw*omegas(1:3);

end

end