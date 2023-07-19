function ds = Dynamics_NP_RW(t, s)
% Note: All angles should be in degrees
% INPUT: s = (theta, thetadot, omegaw, Kp)

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
    
    % Solar Radiation Pressure @ each instant after the eclipse
    Fsrp_s = -PHI/c*(1/Rbepi)^2*A*cosd(saa)*(1-Cs);
    Fsrp_n = -PHI/c*(1/Rbepi)^2*A*cosd(saa)*cosd(saa)*2*Cs;
    Fsrp = sqrt(Fsrp_n^2+Fsrp_s^2);
    delta = acosd(-Fsrp_s/Fsrp);
    Fsrp_ort = Fsrp * cosd(rad2deg(theta) + delta);
    Tsrp = Fsrp_ort * b;
    
    % Gravity Gradient @ each instant after the eclipse
    Rbm0 = (Rm+hca)/cosd(lambda0);
    z0 = -Rbm0*sind(lambda0);
    z = z0 + v0*t;              % position of the SC
    Rbm = sqrt(z^2 + (Rm+hca)^2);
    lambda2 = asind(-z/Rbm);
    Beta2 = rad2deg(theta) + lambda2;
    Gx = 3/2*GMm/(Rbm^3)*(Iz-Iy)*sind(2*Beta2)*cosd(phi0);
    
    Td = Gx+Tsrp;               % value of the total disturbance torque

    % Computation of Desired Pitch Angle and Desired Pitch Angle Derivative
    thetaG = atan(z/(Rm+hca));       
    thetaGdot = -(Rm+hca)*v0/(z^2+(Rm+hca)^2);

    % Computation of Control Torque
    Tc = Kp*(thetaG-theta) + Kd*(thetaGdot-thetadot);

    ds = zeros(length(s), 1);
    
    % Algorythm - NOT Considering Wheel Limitations
    ds(1) = thetadot;
    ds(2) = (-Td + Tc)/Ix;
    ds(3) = -Tc/Iw;    

    
    
    % Algorythm - Considering Wheel Limitations
%     if abs(Tc) < Tmax
% 
%         if abs(omegaw) < wmax
%             ds(1) = thetadot;
%             ds(2) = (-Td+Tc)/Ix;
%             ds(3) = Tc/Iw;
% 
%         elseif (abs(omegaw) >= wmax && (Tc*omegaw) > 0)
%             ds(1) = thetadot;
%             ds(2) = -Td/Ix;
%             ds(3) = 0;
% 
%         elseif (abs(omegaw) >= wmax && (Tc*omegaw) < 0)
%             ds(1) = thetadot;
%             ds(2) = (-Td+Tc)/Ix;
%             ds(3) = Tc/Iw;
% 
%         end
%     end
% 
%     if abs(Tc) >= Tmax
% 
%         if abs(omegaw) < wmax
%             ds(1) = thetadot;
%             ds(2) = (-Td+Tmax*sign(Tc))/Ix;
%             ds(3) = Tmax*sign(Tc)/Iw;
% 
%         elseif (abs(omegaw) >= wmax && (Tc*omegaw) > 0)
%             ds(1) = thetadot;
%             ds(2) = -Td/Ix;
%             ds(3) = 0;
%         
%         elseif (abs(omegaw) >= wmax && (Tc*omegaw) < 0)
%             ds(1) = thetadot;
%             ds(2) = (-Td+Tmax*sign(Tc))/Ix;
%             ds(3) = Tmax*sign(Tc)/Iw;
%         end
%     end

end