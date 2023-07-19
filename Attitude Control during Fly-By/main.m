%% Homework 2 - Leonardo Russo

close all
clear all
clc

%% Data Loading from Known Quantities

A = 42;                 % m^2
PHI = 1371;             % W/m^2
c = 299792458;          % m/s
v0 = 7.7;               % km/s
Rbepi = 0.387;          % AU
saa = 75;               % deg       % sun aspect angle
Cs = 0.3;      
b = 2;                  % m

theta0 = 40.60083;      % deg
psi0 = 0;               % deg
phi0 = 0;               % deg

thetadot0 = 0;          % deg/s
thetaG1 = 40.6;         % deg
w0 = 3000*2*pi()/60;    % rad/s
wmax = 4000*2*pi()/60;  % rad/s

Rm = 2439.7;            % km
hca = 400;              % kme
lambda0 = 80;           % deg
GMm = 22032;            % km^3/s^2

Ix = 8817;              % kg m^2
Iy = 26609;             % kg m^2
Iz = 26971;             % kg m^2

Km = 0.118;             % Nm/A
Kw = 0.003;             % V/rpm
R = 1.2;                % Ohm
Pmax = 29;              % W

Tmax = 0.211;           % Nm

% Utility Vectors
n_errors = 6;
errors = zeros(1, n_errors);
errors_np = zeros(1, n_errors);
fig = 1;

%% Trajectory of the SpaceCraft

% Inertial Pointing Flight
R0 = (Rm+hca)/cosd(lambda0);
z0 = -R0*sind(lambda0);

te1 = -(Rm+z0)/v0;          % time @ which the SC enters the eclipse region
te2 = (Rm-z0)/v0;           % time @ which the SC exits the eclipse region

tmax = 30;                  % s
N = 1000;                   % n° of total time intervals
M = N-tmax-1;               % n° of time intervals after tmax

tspan1 = [0:tmax, linspace((te2+tmax*(M-2))/(M-1), te2, M)]';

z1 = z0 + v0*tspan1;            % position of the SC @ each instant
Rbm1 = sqrt(z1.^2 + (Rm+hca)^2);    % distance between Bepi and Mercury @ each instant
lambda1 = asind(-z1./Rbm1);

% Nadir Pointing
tf = -2*z0/v0;              % final time of the problem for lambda = -80°
tspan_np = linspace(te2, tf, N)';       % same length of tspan1

z2 = z0 + v0*tspan_np;
Rbm2 = sqrt(z2.^2 + (Rm+hca)^2);
lambda2 = asind(-z2./Rbm2);

% Evaluation of Disturbace Torques 
% Note: at every instant I have a different value of disturbance torque

% Solar Radiation Pressure Force
Fsrp_s = -PHI/c*(1/Rbepi)^2*A*cosd(saa)*(1-Cs);
Fsrp_n = -PHI/c*(1/Rbepi)^2*A*cosd(saa)*cosd(saa)*2*Cs;
Fsrp = sqrt(Fsrp_n^2+Fsrp_s^2);
delta = acosd(-Fsrp_s/Fsrp);

%% Attitude Control for Inertial Pointing
% S is the state vector = (theta, thetadot, omega_wheel)

S0 = [deg2rad(theta0), deg2rad(thetadot0), w0];           % state at t = t0

Kp = 160;               % choice of Kp
Kd = 2*sqrt(Kp*Ix);

S = [S0, Kp, deg2rad(thetaG1)]';

Tol0 = 1e-9;
Tol1 = 1e-9;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

eps_max = 1/3600;       % deg

[t, sol] = ode45('Dynamics', tspan1, S, options);

theta = rad2deg(sol(:, 1));
thetadot = rad2deg(sol(:, 2));
w = sol(:, 3);

% Evaluation of Tsrp
Fsrp_ort1 = Fsrp * cosd(theta + delta);
Tsrp1 = Fsrp_ort1 * b;

% Evaluation of Gx
Beta1 = theta+lambda1;
for i = 1:N             % while the SC is in eclipse we will consider Tsrp as negligible
    Gx1(i) = 3/2*GMm/(Rbm1(i)^3)*(Iz-Iy)*sind(2*Beta1(i))*cosd(phi0);   % note that phi is const.
    if tspan1(i) < te1
        Td1(i) = Gx1(i)+Tsrp1(i);           % values of the disturbance torque in each interval
    else
        Td1(i) = Gx1(i);
    end
end

% Evaluation of Wheel Quantities
Tc = Kp*deg2rad(thetaG1-theta) - Kd*deg2rad(thetadot);
iw = Tc/Km;
w_rpm = w*60/(2*pi());
Pw = (Kw*w_rpm + R*iw).*iw;

% Checking Required Conditions

eps = thetaG1 - theta;
over_max_theta = 0;
under_theta_flag = 0;
over_wmax_flag = 0;
over_Pmax_flag = 0;
over_Tmax_flag = 0;

for i = 1:length(tspan1)
    if (abs(eps(i)) < eps_max && under_theta_flag == 0)
        under_theta_flag = 1;
    end
    if (abs(eps(i)) > eps_max && under_theta_flag == 1)
        over_max_theta = over_max_theta + 1;
    end
    if w(i) >= wmax
        over_wmax_flag = 1;
        time_exc_wmax = tspan1(i);
    end

    if Pw(i) >= Pmax
        over_Pmax_flag = 1;
    end

    if abs(Tc(i)) >= Tmax
        over_Tmax_flag = 1;
    end
end

% Pointing Error
if over_max_theta > 0
    fprintf("Error 1: Attitude error does NOT stay in the limits!!!")
    errors(1) = 1;
end

% Maximum time for the Correction
if abs(eps(tspan1(tspan1 == tmax))) >= eps_max
    fprintf("Error 2: Attitude is NOT reached in time!!!")
    errors(2) = 1;
end

% Maximum Wheel Angular Velocity
if over_wmax_flag == 1
    fprintf("Error 3: Exceeded Reaction Wheel Angular Velocity!!!")
    errors(3) = 1;
end

% Maximum Control Torque
if over_Tmax_flag == 1
    fprintf("Error 4: Exceeded Maximum Control Torque!!!")
    errors(4) = 1;
end

% Maximum Power
if over_Pmax_flag == 1
    fprintf("Error 5: Exceeded Maximum Wheel Power!!!")
    errors(5) = 1;
end


% Plots of the Results of the Inertial Pointing

zoom_time = 300;        % s

figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan1, theta, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan1, thetaG1+eps_max*ones(length(tspan1), 1), '--', 'Color', '#2574db', 'LineWidth', 0.5)
plot(tspan1, thetaG1-eps_max*ones(length(tspan1), 1), '--', 'Color', '#2574dc', 'LineWidth', 0.5)
title('Evolution of Pitch Angle')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\theta$ $[deg]$','interpreter','latex','FontSize',14)
% saveas(gcf,'1.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan1, theta, 'Color', '#e62e5f', 'LineWidth', 1.5)
axis([0, zoom_time, thetaG1-eps_max-3*std(theta), theta0+3*std(theta)])
plot(tspan1, thetaG1+eps_max*ones(length(tspan1), 1), '--', 'Color', '#2574db', 'LineWidth', 0.5)
plot(tspan1, thetaG1-eps_max*ones(length(tspan1), 1), '--', 'Color', '#2574dc', 'LineWidth', 0.5)
title('Evolution of Pitch Angle')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\theta$ $[deg]$','interpreter','latex','FontSize',14)
% saveas(gcf,'2.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan1, w, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan1, ones(length(tspan1), 1)*wmax, '--','Color', '#2574db', 'LineWidth', 0.5)
% axis([0, zoom_time, -50, wmax+10*std(w)])
title('Evolution of Wheel Angular Speed')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\omega_{w}$ $[rad/s]$','interpreter','latex','FontSize',14)
% saveas(gcf,'3.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan1, Tc, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan1, ones(length(tspan1), 1)*Tmax, '--','Color', '#2574db', 'LineWidth', 0.5)
% axis([0, zoom_time, -10, Pmax+10*std(Pw)])
title('Evolution of Control Torque')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$T_{c}$ $[Nm]$','interpreter','latex','FontSize',14)
% saveas(gcf,'4.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan1, Pw, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan1, ones(length(tspan1), 1)*Pmax, '--','Color', '#2574db', 'LineWidth', 0.5)
% axis([0, zoom_time, -10, Pmax+10*std(Pw)])
title('Evolution of Wheel Electrical Power')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$P_{w}$ $[W]$','interpreter','latex','FontSize',14)
% saveas(gcf,'5.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan1, iw, 'Color', '#e62e5f', 'LineWidth', 1.5)
% axis([0, zoom_time, -10, Pmax+10*std(Pw)])
title('Evolution of Wheel Current')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$i_{w}$ $[A]$','interpreter','latex','FontSize',14)
% saveas(gcf,'6.jpg')

%% Nadir Pointing with RW

eps_nadir_max = 1;          % deg

Snp0 = [deg2rad(theta(end)), deg2rad(thetadot(end)), w(end)];           % state at t = te2
Snp = [Snp0, Kp]';      % we no longer provide thetaG since it changes @ each instant

[t, sol] = ode45('Dynamics_NP_RW1', tspan_np, Snp, options);

theta_np = rad2deg(sol(:, 1));
thetadot_np = rad2deg(sol(:, 2));
w_np = sol(:, 3);

thetaG_np = atand(z2/(Rm+hca));
thetaGdot_np = -rad2deg((Rm+hca)*v0./(z2.^2+(Rm+hca)^2));

% Evaluation of Tsrp
Fsrp_ort2 = Fsrp * cosd(theta_np + delta);
Tsrp2 = Fsrp_ort2 * b;

% Evaluation of Gx
Beta2 = theta_np + lambda2;
for i = 1:N
    Gx2(i) = 3/2*GMm/(Rbm2(i)^3)*(Iz-Iy)*sind(2*Beta2(i))*cosd(phi0);
end

Td2 = Tsrp2 + Gx2;

Gx = [Gx1; Gx2];
Tsrp = [Tsrp1; Tsrp2];
Td = [Td1; Td2];

% Evaluation of Wheel Quantities
Tc_np = Kp*deg2rad(thetaG_np-theta_np) + Kd*deg2rad(thetaGdot_np-thetadot_np);
iw_np = Tc_np/Km;
w_rpm_np = w_np*60/(2*pi());
Pw_np = (Kw*w_rpm_np + R*iw_np).*iw_np;

% %% Checking all Conditions
% 
% % Check that |theta| < theta_max
% eps_np = thetaG_np - theta_np;
% over_max_theta_np = 0;
% under_theta_flag_np = 0;
% over_wmax_flag_np = 0;
% over_Pmax_flag_np = 0;
% over_Tmax_flag_np = 0;
% 
% for i = 1:length(tspan_np)
%     if (abs(eps_np(i)) < eps_nadir_max && under_theta_flag_np == 0)
%         under_theta_flag_np = 1;
%     end
%     if (abs(eps_np(i)) > eps_nadir_max && under_theta_flag_np == 1)
%         over_max_theta_np = over_max_theta_np + 1;
%     end
%     if w_np(i) >= wmax
%         over_wmax_flag_np = 1;
%         time_exc_wmax_np = tspan_np(i);
%     end
% 
%     if Pw_np(i) >= Pmax
%         over_Pmax_flag_np = 1;
%     end
% 
%     if abs(Tc_np(i)) >= Tmax
%         over_Tmax_flag_np = 1;
%     end
% end
% 
% % Pointing Error
% if over_max_theta_np > 0
%     fprintf("Error 1: Attitude error does NOT stay in the limits during Nadir Pointing with Wheels!!!\n")
%     errors_np(1) = 1;
% end
% 
% % Maximum Wheel Angular Velocity
% if over_wmax_flag_np == 1
%     fprintf("Error 3: Exceeded Reaction Wheel Angular Velocity during Nadir Pointing with Wheels!!!\n")
%     errors_np(3) = 1;
% end
% 
% % Maximum Control Torque
% if over_Tmax_flag_np == 1
%     fprintf("Error 4: Exceeded Maximum Control Torque during Nadir Pointing!!!\n")
%     errors_np(4) = 1;
% end
% 
% % Maximum Power
% if over_Pmax_flag_np == 1
%     fprintf("Error 5: Exceeded Maximum Wheel Power during Nadir Pointing!!!\n")
%     errors_np(5) = 1;
% end


%% Plots of the Results of Nadir Pointing using RW

zoom_time_np = te2+zoom_time;        % s

figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan_np, theta_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan_np, thetaG_np+eps_nadir_max*ones(length(tspan_np), 1)', '--', 'Color', '#2574db', 'LineWidth', 0.5)
plot(tspan_np, thetaG_np-eps_nadir_max*ones(length(tspan_np), 1)', '--', 'Color', '#2574dc', 'LineWidth', 0.5)
title('Evolution of Pitch Angle')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\theta$ $[deg]$','interpreter','latex','FontSize',14)
% saveas(gcf,'7.jpg')


% figure(fig)
% fig = fig+1;
% 
% set(gcf, 'Position',  [850, 150, 650, 400])
% hold on
% plot(tspan_np, theta_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
% axis([te2, zoom_time_np, theta_np(1)-std(theta_np), max(thetaG_np(tspan_np<zoom_time_np))])
% plot(tspan_np, thetaG_np+eps_nadir_max*ones(length(tspan_np), 1)', '--', 'Color', '#2574db', 'LineWidth', 0.5)
% plot(tspan_np, thetaG_np-eps_nadir_max*ones(length(tspan_np), 1)', '--', 'Color', '#2574dc', 'LineWidth', 0.5)
% title('Evolution of Pitch Angle')
% xlabel('t\,[s]','Interpreter','latex','FontSize',12)
% ylabel('$\theta$ $[deg]$','interpreter','latex','FontSize',14)
% 
% 
% figure(fig)
% fig = fig+1;
% 
% set(gcf, 'Position',  [850, 150, 650, 400])
% hold on
% plot(tspan_np, w_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
% plot(tspan_np, wmax*ones(length(tspan_np), 1)', '--','Color', '#2574db', 'LineWidth', 0.5)
% % axis([te2, te2+10, w_np(1)-10, wmax+3*std(w_np)]) this is good as it is
% title('Evolution of Wheel Angular Speed')
% xlabel('t\,[s]','Interpreter','latex','FontSize',12)
% ylabel('$\omega_{w}$ $[rad/s]$','interpreter','latex','FontSize',14)
% 
% 
% figure(fig)
% fig = fig+1;
% 
% set(gcf, 'Position',  [850, 150, 650, 400])
% hold on
% plot(tspan_np, Tc_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
% plot(tspan_np, Tmax*ones(length(tspan_np), 1)', '--','Color', '#2574db', 'LineWidth', 0.5)
% % axis([te2, te2+10, -Tmax, max(Tc_np(tspan_np<te2+10))])
% title('Evolution of Control Torque')
% xlabel('t\,[s]','Interpreter','latex','FontSize',12)
% ylabel('$T_{c}$ $[Nm]$','interpreter','latex','FontSize',14)
% 
% 
% figure(fig)
% fig = fig+1;
% 
% set(gcf, 'Position',  [850, 150, 650, 400])
% hold on
% plot(tspan_np, Pw_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
% plot(tspan_np, Pmax*ones(length(tspan_np), 1)', '--','Color', '#2574db', 'LineWidth', 0.5)
% % axis([te2, zoom_time_np, -10, Pmax+10*std(Pw)])
% title('Evolution of Wheel Electrical Power')
% xlabel('t\,[s]','Interpreter','latex','FontSize',12)
% ylabel('$P_{w}$ $[W]$','interpreter','latex','FontSize',14)

%% Nadir Pointing with RW + Thrusters

eps_nadir_max = 1;          % deg

Snp0 = [deg2rad(theta(end)), deg2rad(thetadot(end)), w(end)];           % state at t = te2
Snp = [Snp0, Kp]';      % we no longer provide thetaG since it changes @ each instant

tstart = te2;
t_thrust = 0.7;             % s
tend = te2 + t_thrust;

[t1, sol1] = ode45('Dynamics_NP_RWT', [tstart, tend], Snp, options);

theta_np = rad2deg(sol1(:, 1));
thetadot_np = rad2deg(sol1(:, 2));
w_np = sol1(:, 3);

Snp = [sol1(end, 1), sol1(end, 2), sol1(end, 3), Kp]';

tstart = tend;
tend = tf;

[t2, sol] = ode45('Dynamics_NP_RW', linspace(tstart, tend, N), Snp, options);

theta_np = [theta_np; rad2deg(sol(:, 1))];
thetadot_np = [thetadot_np; rad2deg(sol(:, 2))];
w_np = [w_np; sol(:, 3)];

tspan2 = [t1; t2];

z2 = z0 + v0*tspan2;
Rbm2 = sqrt(z2.^2 + (Rm+hca)^2);
lambda2 = asind(-z2./Rbm2);

thetaG_np = atand(z2/(Rm+hca));
thetaGdot_np = -rad2deg((Rm+hca)*v0./(z2.^2+(Rm+hca)^2));

% % Evaluation of Tsrp
% Fsrp_ort2 = Fsrp * cosd(theta_np + delta);
% Tsrp2 = Fsrp_ort2 * b;
% 
% % Evaluation of Gx
% Beta2 = theta_np + lambda2;
% for i = 1:N
%     Gx2(i) = 3/2*GMm/(Rbm2(i)^3)*(Iz-Iy)*sind(2*Beta2(i))*cosd(phi0);
% end
% 
% Td2 = Tsrp2 + Gx2;
% 
% Gx = [Gx1; Gx2];
% Tsrp = [Tsrp1; Tsrp2];
% Td = [Td1; Td2];

% Evaluation of Wheel Quantities
Tc_np = Kp*deg2rad(thetaG_np-theta_np) + Kd*deg2rad(thetaGdot_np-thetadot_np);
for i=1:length(Tc_np)
    if abs(Tc_np(i)) > Tmax
        Tc_np(i) = Tmax*sign(Tc_np(i));
    end
    if abs(w_np(i)) > wmax
        w_np(i) = wmax*sign(w_np(i));
    end
end
iw_np = Tc_np/Km;
w_rpm_np = w_np*60/(2*pi());
Pw_np = (Kw*w_rpm_np + R*iw_np).*iw_np;

%% Checking all Conditions
% 
% % Check that |theta| < theta_max
% eps_np = thetaG_np - theta_np;
% over_max_theta_np = 0;
% under_theta_flag_np = 0;
% over_wmax_flag_np = 0;
% over_Pmax_flag_np = 0;
% over_Tmax_flag_np = 0;
% 
% for i = 1:length(tspan2)
%     if (abs(eps_np(i)) < eps_nadir_max && under_theta_flag_np == 0)
%         under_theta_flag_np = 1;
%     end
%     if (abs(eps_np(i)) > eps_nadir_max && under_theta_flag_np == 1)
%         over_max_theta_np = over_max_theta_np + 1;
%     end
%     if w_np(i) >= wmax
%         over_wmax_flag_np = 1;
%         time_exc_wmax_np = tspan2(i);
%     end
% 
%     if Pw_np(i) >= Pmax
%         over_Pmax_flag_np = 1;
%     end
% 
%     if abs(Tc_np(i)) >= Tmax
%         over_Tmax_flag_np = 1;
%     end
% end
% 
% % Pointing Error
% if over_max_theta_np > 0
%     fprintf("Error 1: Attitude error does NOT stay in the limits during Nadir Pointing!!!\n")
%     errors_np(1) = 1;
% end
% 
% % Maximum Wheel Angular Velocity
% if over_wmax_flag_np == 1
%     fprintf("Error 3: Exceeded Reaction Wheel Angular Velocity during Nadir Pointing!!!\n")
%     errors_np(3) = 1;
% end
% 
% % Maximum Control Torque
% if over_Tmax_flag_np == 1
%     fprintf("Error 4: Exceeded Maximum Control Torque during Nadir Pointing!!!\n")
%     errors_np(4) = 1;
% end
% 
% % Maximum Power
% if over_Pmax_flag_np == 1
%     fprintf("Error 5: Exceeded Maximum Wheel Power during Nadir Pointing!!!\n")
%     errors_np(5) = 1;
% end


%% Plots of the Results of Nadir Pointing using RW + Thrusters

zoom_time = 300; 

zoom_time_np = te2+zoom_time;        % s

figure(fig)
fig = fig+1;

thetaG_up = thetaG_np+eps_nadir_max*ones(length(tspan2), 1)';
thetaG_down = thetaG_np-eps_nadir_max*ones(length(tspan2), 1)';

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan2, theta_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan2, thetaG_np+eps_nadir_max*ones(length(tspan2), 1)', '--', 'Color', '#2574db', 'LineWidth', 0.5)
plot(tspan2, thetaG_np-eps_nadir_max*ones(length(tspan2), 1)', '--', 'Color', '#2574dc', 'LineWidth', 0.5)
title('Evolution of Pitch Angle')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\theta$ $[deg]$','interpreter','latex','FontSize',14)
% saveas(gcf,'8.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan2, theta_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
axis([te2, zoom_time_np, theta_np(1)-std(theta_np), max(thetaG_np(tspan2<zoom_time_np))])
plot(tspan2, thetaG_np+eps_nadir_max*ones(length(tspan2), 1)', '--', 'Color', '#2574db', 'LineWidth', 0.5)
plot(tspan2, thetaG_np-eps_nadir_max*ones(length(tspan2), 1)', '--', 'Color', '#2574dc', 'LineWidth', 0.5)
title('Evolution of Pitch Angle')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\theta$ $[deg]$','interpreter','latex','FontSize',14)
% saveas(gcf,'9.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan2, w_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan2, wmax*ones(length(tspan2), 1)', '--','Color', '#2574db', 'LineWidth', 0.5)
% axis([te2, te2+10, w_np(1)-10, wmax+3*std(w_np)]) this is good as it is
title('Evolution of Wheel Angular Speed')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\omega_{w}$ $[rad/s]$','interpreter','latex','FontSize',14)
% saveas(gcf,'10.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan2, Tc_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan2, Tmax*ones(length(tspan2), 1)', '--','Color', '#2574db', 'LineWidth', 0.5)
% axis([te2, te2+10, -Tmax, max(Tc_np(tspan2<te2+10))])
title('Evolution of Control Torque')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$T_{c}$ $[Nm]$','interpreter','latex','FontSize',14)
% saveas(gcf,'11.jpg')


figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan2, Pw_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
plot(tspan2, Pmax*ones(length(tspan2), 1)', '--','Color', '#2574db', 'LineWidth', 0.5)
% axis([te2, zoom_time_np, -10, Pmax+10*std(Pw)])
title('Evolution of Wheel Electrical Power')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$P_{w}$ $[W]$','interpreter','latex','FontSize',14)
% saveas(gcf,'12.jpg')

figure(fig)
fig = fig+1;

set(gcf, 'Position',  [850, 150, 650, 400])
hold on
plot(tspan2, iw_np, 'Color', '#e62e5f', 'LineWidth', 1.5)
% axis([te2, zoom_time_np, -10, Pmax+10*std(Pw)])
title('Evolution of Wheel Current')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$i_{w}$ $[A]$','interpreter','latex','FontSize',14)
% saveas(gcf,'13.jpg')

%% Desaturation Manoeuver

thetaGn = 51.5;         % deg




