function [R] = RotECI2Peri(COE)
% Description: this function generates a rotation matrix from the Inertial
% Frame to the (r, theta, h) frame from the COE given as input.
% https://en.wikiversity.org/wiki/PlanetPhysics/Euler_313_Sequence

incl = COE(3);
Omega = COE(4);
omega = COE(5);
theta_star = COE(6);

theta_t = omega + theta_star;

R3_Omega = [cos(Omega) sin(Omega) 0;...
            -sin(Omega) cos(Omega) 0;...
            0 0 1];

R1_incl = [1 0 0;...
           0 cos(incl) sin(incl);...
           0 -sin(incl) cos(incl)];

R3_theta_t = [cos(theta_t) sin(theta_t) 0;...
              -sin(theta_t) cos(theta_t) 0;...
              0 0 1];

R = R3_theta_t * R1_incl * R3_Omega;

end