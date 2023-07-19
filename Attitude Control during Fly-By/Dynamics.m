function ds = Dynamics(t, s)
% Input: s = (theta, thetadot, omegaw, Kp,  thetaG)

    % Report Known Quantities
    A = 42;             % m^2
    PHI = 1371;         % W/m^2
    c = 299792458;      % m/s
    v0 = 7.7;           % km/s
    Rbepi = 0.387;      % AU
    saa = 75;           % deg       % sun aspect angle
    Cs = 0.3;      
    b = 2;              % m
    Rm = 2439.7;        % km
    hca = 400;          % km
    lambda0 = 80;       % deg
    GMm = 22032;        % km^3/s^2
    Ix = 8817;          % kg m^2
    Iy = 26609;         % kg m^2
    Iz = 26971;         % kg m^2
    Iw = 0.054;         % kg m^2
    phi0 = 0;           % deg
    Tmax = 0.211;           % Nm
    wmax = 4000*2*pi()/60;  % rad/s
    Tth = 40;           % Nm

    % Reading input values
    theta = s(1);       % rad
    thetadot = s(2);    % rad
    omegaw = s(3);      % rad/s
    Kp = s(4);
    Kd = 2*sqrt(Kp*Ix);
    thetaG = s(5);      % rad

    % Solar Radiation Pressure
    Fsrp_s = -PHI/c*(1/Rbepi)^2*A*cosd(saa)*(1-Cs);
    Fsrp_n = -PHI/c*(1/Rbepi)^2*A*cosd(saa)*cosd(saa)*2*Cs;
    Fsrp = sqrt(Fsrp_n^2+Fsrp_s^2);
    delta = acosd(-Fsrp_s/Fsrp);
    Fsrp_ort = Fsrp * cosd(rad2deg(theta) + delta);
    Tsrp = Fsrp_ort * b;

    % Gravity Gradient
    Rbm0 = (Rm+hca)/cosd(lambda0);
    z0 = -Rbm0*sind(lambda0);
    z = z0 + v0*t;          % position of the SC @ every instant
    Rbm = sqrt(z^2 + (Rm+hca)^2);
    lambda1 = asind(-z/Rbm);
    Beta1 = rad2deg(theta) + lambda1;
    Gx = 3/2*GMm/(Rbm^3)*(Iz-Iy)*sind(2*Beta1)*cosd(phi0);
    
    te1 = -(Rm+z0)/v0;          % time @ which the SC enters the eclipse region
    te2 = (Rm-z0)/v0;           % time @ which the SC exits the eclipse region

    if t <= te1
        Td = Gx+Tsrp;           % values of the disturbance torque at each instant before eclipse
    elseif t >te1 && t <= te2
        Td = Gx;
    end    

    ds = zeros(length(s), 1);

    Tc = Kp*(thetaG-theta) - Kd*thetadot;

    % Td is positive but it acts opposite to omegaw

%     % Algorythm - NOT Considering Wheel Limitations
%     ds(1) = thetadot;
%     ds(2) = (Td + Tc)/Ix;
%     ds(3) = -Tc/Iw;

    % Algorythm - Considering Wheel Limitations
    if abs(Tc) < Tmax

        if abs(omegaw) < wmax
            ds(1) = thetadot;
            ds(2) = (-Td+Tc)/Ix;
            ds(3) = Tc/Iw;

        elseif (abs(omegaw) >= wmax && (Tc*omegaw) < 0)
            ds(1) = thetadot;
            ds(2) = -Td/Ix;
            ds(3) = 0;

        elseif (abs(omegaw) >= wmax && (Tc*omegaw) > 0)
            ds(1) = thetadot;
            ds(2) = (-Td+Tc)/Ix;
            ds(3) = Tc/Iw;

        end
    end

    if abs(Tc) >= Tmax

        if abs(omegaw) < wmax
            ds(1) = thetadot;
            ds(2) = (-Td+Tmax*sign(Tc))/Ix;
            ds(3) = Tmax*sign(Tc)/Iw;

        elseif (abs(omegaw) >= wmax && (Tc*omegaw) < 0)
            ds(1) = thetadot;
            ds(2) = -Td/Ix;
            ds(3) = 0;
        
        elseif (abs(omegaw) >= wmax && (Tc*omegaw) > 0)
            ds(1) = thetadot;
            ds(2) = (-Td+Tmax*sign(Tc))/Ix;
            ds(3) = Tmax*sign(Tc)/Iw;
        end
    end
end