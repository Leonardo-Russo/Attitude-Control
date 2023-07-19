%% Homework 3 - Leonardo Russo 2015563

close all
clear all
clc

addpath('Library/')

%% Options

% Define the options for the Integration
Tol0 = 1e-11;
Tol1 = 1e-13;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

savechoice = '1';

optional_plots = '1';


%% Define Known Quantities

global mu J Is sign_qe0_0 section omega_dot_max

mu = 398600.4415;       % km^3/s^2
D_sid = 86164;          % s
T_Mol = D_sid / 2;      % s

omega_dot_max = deg2rad(30);       % rad/s^2

Nw = 4;     % n° of RW

Jc = [900, 50, -100;...
      50, 1100, 150;...
      -100, 150, 1250];     % kg m^2

Is = 1;     % kg m^2
It = 0.5;   % kg m^2

J = Jc;     % assuming J_Wheels << J_S/C

omega1_0 = 0.5;     % rad/s
omega2_0 = 0.5;     % rad/s
omega3_0 = -0.5;    % rad/s
omega4_0 = -0.5;    % rad/s
omegas_0 = [omega1_0, omega2_0, omega3_0, omega4_0]';

qb_0 = [0.1, 0.3, -0.5]';
qb0_0 = sqrt(1 - norm(qb_0)^2);
xb0 = [qb0_0; qb_0];

w_0 = [-0.1, 0.05, 0]';     % rad/s


% Orbital Parameters
a = (mu * (T_Mol/(2*pi))^2)^(1/3);    % km
e = 0.73;
incl = 63.4;        % deg
Omega = 0;          % deg
w_i = -90;          % deg
theta_star = 90;    % deg        
COE_0 = [a, e, deg2rad(incl), deg2rad(Omega), deg2rad(w_i), deg2rad(theta_star)];

[r0_vect, v0_vect] = coe2rvECI(COE_0, mu);


% Define Commanded Attitude
R_B = [0 1 0;...
       0 0 -1;...
       -1 0 0];     % Commanded Attitude from Perifocal f.o.r.

% Define Rotation Matrix from ECI to Perifocal f.o.r.
RNP = RotECI2Peri(COE_0);

[qc0_0, qc_0] = C2q(R_B*RNP);

% Define Error Quaternions
qe0_0 = qc0_0*qb0_0 + qc_0' * qb_0;
qe_0 = -qc_0*qb0_0 + qc0_0*qb_0 - skew(qc_0)*qb_0;
xe0 = [qe0_0; qe_0];

sign_qe0_0 = sign(qe0_0);


% Gain Parameters
omega_n = 0.03;     % rad/s
xi = 1;
c1 = 2 * omega_n^2;
c2 = 2*xi*omega_n/c1;
invA = c1*eye(3);
B = c2 * eye(3);


% Define aj as the axis of each RW
a1 = [1 0 0]';
a2 = [0 1 0]';
a3 = [0 0 1]';
a4 = [sqrt(3)/3 sqrt(3)/3 sqrt(3)/3]';
as = [a1, a2, a3, a4];

% Compute Constant Variables
Mc = 0;
Ja = buildJa(Is, as);


%% Section A

t0 = 0;         % sec
tf = 1200;      % sec
section = 'A';

X0 = [w_0; omegas_0; r0_vect; v0_vect; xe0];

[tspan, X] = ode113('DynamicalModel', [t0, tf], X0, options);

M = length(tspan);


% Retrieve Data from Propagation
w = X(:, 1:3);
omegas = X(:, 4:7);
r_vect = X(:, 8:10);
v_vect = X(:, 11:13);
xe = X(:, 14:17);

qe0 = xe(:, 1);
qe = xe(:, 2:4);


% Compute the Variables inside the Dynamical Model
we = zeros(M, 3);
Tc = zeros(M, 3);
Ta = zeros(M, 3);
omegas_dot = zeros(M, 4);


for i = 1 : M   % repeats the same steps followed inside Dynamical Model

    w_i = w(i, :)';
    omegas_i = omegas(i, :)';
    r_vect_i = r_vect(i, :)';
    v_vect_i = v_vect(i, :)';
    qe0_i = qe0(i);
    qe_i = qe(i, :)';

    r = norm(r_vect_i);
    h = norm(cross(r_vect_i, v_vect_i));

    COE = rvECI2coe(r_vect_i, v_vect_i, mu);
    RNP = RotECI2Peri(COE);         % Rotation Matrix from N to Perifocal

    v_vectP = RNP * v_vect_i;
    vr = v_vectP(1);            % Radial Velocity

    wc = [0 -h/r^2 0]';
    wc_dot = [0 2*h/r^3*vr 0]';
    
    RCB = q2C(qe0_i, qe_i);
    we_i = w_i - RCB * wc;
    wd = w_i - wc;

    Tc_i = skew(w_i)*J*w_i - Mc + J*wc_dot - J*invA*B*wd - sign_qe0_0*J*invA*qe_i;
    Ta_i = Tc_i - skew(w_i)*Ja*omegas_i;
    omegas_dot_i = -Ja' * inv(Ja * Ja') * Tc_i;

    % Assign Computed Values
    we(i, :) = we_i;
    Tc(i, :) = Tc_i;
    Ta(i, :) = Ta_i;
    omegas_dot(i, :) = omegas_dot_i;

end


%%% Plot of the Results of Section A %%%

figure('Name', "Error Quaternions in Section A")

hold on
grid on
plot(tspan, qe0)
plot(tspan, qe(:, 1))
plot(tspan, qe(:, 2))
plot(tspan, we(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{q_{e,i}}$', 'Interpreter','latex', 'FontSize', 12)
legend('$q_{e0}$', '$q_{e,1}$', '$q_{e,2}$', '$q_{e,3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\quatA.jpg'))
end


figure('Name', "Error Angular Rate in Section A")

hold on
grid on
plot(tspan, we(:, 1))
plot(tspan, we(:, 2))
plot(tspan, we(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\omega_{e,i}} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{e,1}$', '$\omega_{e,2}$', '$\omega_{e,3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\weA.jpg'))
end


figure('Name', "Commanded Torque in Section A")

hold on
grid on
plot(tspan, Tc(:, 1))
plot(tspan, Tc(:, 2))
plot(tspan, Tc(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{T_c} \ [Nm]$', 'Interpreter','latex', 'FontSize', 12)
legend('$T_{c1}$', '$T_{c2}$', '$T_{c3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\TcA.jpg'))
end


figure('Name', "Applied Torque in Section A")

hold on
grid on
plot(tspan, Ta(:, 1))
plot(tspan, Ta(:, 2))
plot(tspan, Ta(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{T_a} \ [Nm]$', 'Interpreter','latex', 'FontSize', 12)
legend('$T_{a1}$', '$T_{a2}$', '$T_{a3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\TaA.jpg'))
end


omegas = rad2deg(omegas);
figure('Name', "Wheels Angular Velocity in Section A")

hold on
grid on
plot(tspan, omegas(:, 1))
plot(tspan, omegas(:, 2))
plot(tspan, omegas(:, 3))
plot(tspan, omegas(:, 4))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\omega_{s,i}} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{s1}$', '$\omega_{s2}$', '$\omega_{s3}$', '$\omega_{s4}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\omegasA.jpg'))
end


omegas_dot = rad2deg(omegas_dot);
figure('Name', "Wheels Angular Acceleration in Section A")

hold on
grid on
plot(tspan, omegas_dot(:, 1))
plot(tspan, omegas_dot(:, 2))
plot(tspan, omegas_dot(:, 3))
plot(tspan, omegas_dot(:, 4))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\dot{\omega_{s,i}}} \ [deg/s^2]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\dot{\omega_{s1}}$', '$\dot{\omega_{s2}}$', '$\dot{\omega_{s3}}$', '$\dot{\omega_{s4}}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\omegasdotA.jpg'))
end


if optional_plots == '1'

    figure('Name', "Body Angular Velocity in Section A")
    
    hold on
    plot(tspan, w(:, 1))
    plot(tspan, w(:, 2))
    plot(tspan, w(:, 3))
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('$\underline{\omega} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$', 'Location','best', 'interpreter', 'latex')
    if savechoice == '1'
        saveas(gcf, strcat('Output\omegaA.jpg'))
    end

end


%% Section B

t0 = 0;         % sec
tf = 1 * T_Mol;     % sec
section = 'B';

X0 = [w_0; omegas_0; r0_vect; v0_vect; xe0];

[tspan, X] = ode113('DynamicalModel', [t0, tf], X0, options);

M = length(tspan);


% Retrieve Data from Propagation
w = X(:, 1:3);
omegas = X(:, 4:7);
r_vect = X(:, 8:10);
v_vect = X(:, 11:13);
xe = X(:, 14:17);

qe0 = xe(:, 1);
qe = xe(:, 2:4);


% Compute the Variables inside the Dynamical Model
we = zeros(M, 3);
Tc = zeros(M, 3);
Ta = zeros(M, 3);
omegas_dot = zeros(M, 4);


for i = 1 : M   % repeats the same steps followed inside Dynamical Model

    w_i = w(i, :)';
    omegas_i = omegas(i, :)';
    r_vect_i = r_vect(i, :)';
    v_vect_i = v_vect(i, :)';
    qe0_i = qe0(i);
    qe_i = qe(i, :)';

    r = norm(r_vect_i);
    h = norm(cross(r_vect_i, v_vect_i));

    COE = rvECI2coe(r_vect_i, v_vect_i, mu);
    RNP = RotECI2Peri(COE);         % Rotation Matrix from N to Perifocal

    v_vectP = RNP * v_vect_i;
    vr = v_vectP(1);            % Radial Velocity

    wc = [0 -h/r^2 0]';
    wc_dot = [0 2*h/r^3*vr 0]';
    
    RCB = q2C(qe0_i, qe_i);
    we_i = w_i - RCB * wc;
    wd = w_i - wc;

    Tc_i = skew(w_i)*J*w_i - Mc + J*wc_dot - J*invA*B*wd - sign_qe0_0*J*invA*qe_i;
    Ta_i = Tc_i - skew(w_i)*Ja*omegas_i;
    omegas_dot_i = -Ja' * inv(Ja * Ja') * Tc_i;

    % Assign Computed Values
    we(i, :) = we_i;
    Tc(i, :) = Tc_i;
    Ta(i, :) = Ta_i;
    omegas_dot(i, :) = omegas_dot_i;

end

DrawTraj_3D(r_vect, xe)


%%% Plot of the Results of Section B %%%

figure('Name', "Error Quaternions in Section B")

hold on
grid on
plot(tspan, qe0)
plot(tspan, qe(:, 1))
plot(tspan, qe(:, 2))
plot(tspan, we(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{q_{e,i}}$', 'Interpreter','latex', 'FontSize', 12)
legend('$q_{e0}$', '$q_{e,1}$', '$q_{e,2}$', '$q_{e,3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\quatB.jpg'))
end


figure('Name', "Error Angular Rate in Section B")

hold on
grid on
plot(tspan, we(:, 1))
plot(tspan, we(:, 2))
plot(tspan, we(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\omega_{e,i}} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{e,1}$', '$\omega_{e,2}$', '$\omega_{e,3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\weB.jpg'))
end


figure('Name', "Commanded Torque in Section B")

hold on
grid on
plot(tspan, Tc(:, 1))
plot(tspan, Tc(:, 2))
plot(tspan, Tc(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{T_c} \ [Nm]$', 'Interpreter','latex', 'FontSize', 12)
legend('$T_{c1}$', '$T_{c2}$', '$T_{c3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\TcB.jpg'))
end


if optional_plots == '1'

    figure('Name', "Applied Torque in Section B")

    hold on
    grid on
    plot(tspan, Ta(:, 1))
    plot(tspan, Ta(:, 2))
    plot(tspan, Ta(:, 3))
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('$\underline{T_a} \ [Nm]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$T_{a1}$', '$T_{a2}$', '$T_{a3}$', 'Location','best', 'interpreter', 'latex')
    if savechoice == '1'
        saveas(gcf, strcat('Output\TaB.jpg'))
    end
    
    
    omegas = rad2deg(omegas);
    figure('Name', "Wheels Angular Velocity in Section B")
    
    hold on
    grid on
    plot(tspan, omegas(:, 1))
    plot(tspan, omegas(:, 2))
    plot(tspan, omegas(:, 3))
    plot(tspan, omegas(:, 4))
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('$\underline{\omega_{s,i}} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$\omega_{s1}$', '$\omega_{s2}$', '$\omega_{s3}$', '$\omega_{s4}$', 'Location','best', 'interpreter', 'latex')
    if savechoice == '1'
        saveas(gcf, strcat('Output\omegasB.jpg'))
    end
    
    
    omegas_dot = rad2deg(omegas_dot);
    figure('Name', "Wheels Angular Acceleration in Section B")
    
    hold on
    grid on
    plot(tspan, omegas_dot(:, 1))
    plot(tspan, omegas_dot(:, 2))
    plot(tspan, omegas_dot(:, 3))
    plot(tspan, omegas_dot(:, 4))
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('$\underline{\dot{\omega_{s,i}}} \ [deg/s^2]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$\dot{\omega_{s1}}$', '$\dot{\omega_{s2}}$', '$\dot{\omega_{s3}}$', '$\dot{\omega_{s4}}$', 'Location','best', 'interpreter', 'latex')
    if savechoice == '1'
        saveas(gcf, strcat('Output\omegasdotB.jpg'))
    end


    figure('Name', "Body Angular Velocity in Section B")
    
    hold on
    plot(tspan, w(:, 1))
    plot(tspan, w(:, 2))
    plot(tspan, w(:, 3))
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('$\underline{\omega} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$', 'Location','best', 'interpreter', 'latex')
    if savechoice == '1'
        saveas(gcf, strcat('Output\omegaB.jpg'))
    end

end


%% Section C

t0 = 0;         % sec
tf = 1200;      % sec
section = 'C';

X0 = [w_0; omegas_0; r0_vect; v0_vect; xe0];

[tspan, X] = ode113('DynamicalModel', [t0, tf], X0, options);

M = length(tspan);


% Retrieve Data from Propagation
w = X(:, 1:3);
omegas = X(:, 4:7);
r_vect = X(:, 8:10);
v_vect = X(:, 11:13);
xe = X(:, 14:17);

qe0 = xe(:, 1);
qe = xe(:, 2:4);


% Compute the Variables inside the Dynamical Model
we = zeros(M, 3);
Tc = zeros(M, 3);
Ta = zeros(M, 3);
omegas_dot = zeros(M, 4);


for i = 1 : M   % repeats the same steps followed inside Dynamical Model

    w_i = w(i, :)';
    omegas_i = omegas(i, :)';
    r_vect_i = r_vect(i, :)';
    v_vect_i = v_vect(i, :)';
    qe0_i = qe0(i);
    qe_i = qe(i, :)';

    r = norm(r_vect_i);
    h = norm(cross(r_vect_i, v_vect_i));

    COE = rvECI2coe(r_vect_i, v_vect_i, mu);
    RNP = RotECI2Peri(COE);         % Rotation Matrix from N to Perifocal

    v_vectP = RNP * v_vect_i;
    vr = v_vectP(1);            % Radial Velocity

    wc = [0 -h/r^2 0]';
    wc_dot = [0 2*h/r^3*vr 0]';
    
    RCB = q2C(qe0_i, qe_i);
    we_i = w_i - RCB * wc;
    wd = w_i - wc;

    Tc_i = skew(w_i)*J*w_i - Mc + J*wc_dot - J*invA*B*wd - sign_qe0_0*J*invA*qe_i;
    Ta_i = Tc_i - skew(w_i)*Ja*omegas_i;
    upd = -Ja' * inv(Ja * Ja') * Tc_i;
    
    omegas_dot_i = zeros(4, 1);
    for j = 1 : 4
        if upd(j) > omega_dot_max || upd(j) < -omega_dot_max
            omegas_dot_i(j) = omega_dot_max * sign(upd(j));
        else
            omegas_dot_i(j) = upd(j);
        end
    end

    % Assign Computed Values
    we(i, :) = we_i;
    Tc(i, :) = Tc_i;
    Ta(i, :) = Ta_i;
    omegas_dot(i, :) = omegas_dot_i;

end


%%% Plot of the Results of Section C %%%

figure('Name', "Error Quaternions in Section C")

hold on
grid on
plot(tspan, qe0)
plot(tspan, qe(:, 1))
plot(tspan, qe(:, 2))
plot(tspan, we(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{q_{e,i}}$', 'Interpreter','latex', 'FontSize', 12)
legend('$q_{e0}$', '$q_{e,1}$', '$q_{e,2}$', '$q_{e,3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\quatC.jpg'))
end


figure('Name', "Error Angular Rate in Section C")

hold on
grid on
plot(tspan, we(:, 1))
plot(tspan, we(:, 2))
plot(tspan, we(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\omega_{e,i}} \ [rad/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{e,1}$', '$\omega_{e,2}$', '$\omega_{e,3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\weC.jpg'))
end


figure('Name', "Commanded Torque in Section C")

hold on
grid on
plot(tspan, Tc(:, 1))
plot(tspan, Tc(:, 2))
plot(tspan, Tc(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{T_c} \ [Nm]$', 'Interpreter','latex', 'FontSize', 12)
legend('$T_{c1}$', '$T_{c2}$', '$T_{c3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\TcC.jpg'))
end


figure('Name', "Applied Torque in Section C")

hold on
grid on
plot(tspan, Ta(:, 1))
plot(tspan, Ta(:, 2))
plot(tspan, Ta(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{T_a} \ [Nm]$', 'Interpreter','latex', 'FontSize', 12)
legend('$T_{a1}$', '$T_{a2}$', '$T_{a3}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\TaC.jpg'))
end


omegas = rad2deg(omegas);
figure('Name', "Wheels Angular Velocity in Section C")

hold on
grid on
plot(tspan, omegas(:, 1))
plot(tspan, omegas(:, 2))
plot(tspan, omegas(:, 3))
plot(tspan, omegas(:, 4))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\omega_{s,i}} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\omega_{s1}$', '$\omega_{s2}$', '$\omega_{s3}$', '$\omega_{s4}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\omegasC.jpg'))
end


omegas_dot = rad2deg(omegas_dot);
figure('Name', "Wheels Angular Acceleration in Section C")

hold on
grid on
plot(tspan, omegas_dot(:, 1))
plot(tspan, omegas_dot(:, 2))
plot(tspan, omegas_dot(:, 3))
plot(tspan, omegas_dot(:, 4))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\dot{\omega_{s,i}}} \ [deg/s^2]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\dot{\omega_{s1}}$', '$\dot{\omega_{s2}}$', '$\dot{\omega_{s3}}$', '$\dot{\omega_{s4}}$', 'Location','best', 'interpreter', 'latex')
if savechoice == '1'
    saveas(gcf, strcat('Output\omegasdotC.jpg'))
end


if optional_plots == '1'

    figure('Name', "Body Angular Velocity in Section C")
    
    hold on
    plot(tspan, w(:, 1))
    plot(tspan, w(:, 2))
    plot(tspan, w(:, 3))
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('$\underline{\omega} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$', 'Location','best', 'interpreter', 'latex')
    if savechoice == '1'
        saveas(gcf, strcat('Output\omegaC.jpg'))
    end

end




