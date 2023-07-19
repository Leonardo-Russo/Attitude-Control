%% Homework 1 - Leonardo Russo 2015563

% Note: All spatial quantities will be represented in the Body Frame

close all
clear all
clc

addpath("Library\")

%% Define Known Quantities

% Please, select one of the two cases
global choice

choice = char(string(input('Please, select the case:\n1. Case 1\n2. Case 2\n')));
clc

savechoice = char(string(input('Do you wish to save a copy of the graphs locally?\n1. Yes\n2. No\n')));
clc

mb = 100;       % kg
md = 10;        % kg

% Assume three equal Wheels
mw = 5;         % kg
Rw = 1;         % m
bw = 2;         % m

M = mb + 3*mw + md;     % total mass M

Is = mw * Rw^2 / 2;     % kg m^2
It = mw * Rw^2 / 4;     % kg m^2

Kd = 5.5;       % kg/s^2
Cd = 30;        % kg/s
Kw = 0.1;       % kg m^2/s

% Define Time Domain
t0 = 0;         % s
tf = 30 * 60;   % s

%% Initial Conditions

vp0 = zeros(3, 1);              % m/s
omega0 = deg2rad([36 3 3]');    % rad/s
ws10 = deg2rad(720);            % rad/s
ws20 = deg2rad(60);             % rad/s
ws30 = deg2rad(720);            % rad/s
xi0 = 0;
xidot0 = 0;

nu0 = [vp0; omega0; ws10; ws20; ws30; xidot0];     % IC for nu

y0 = [nu0; xi0];     % IC for y

%% Propagation

% Define the Options for the Integration
Tol0 = 1e-11;
Tol1 = 1e-13;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

[tspan, y] = ode113('Dynamical_Model', [t0, tf], y0, options);

nu = y(:, 1:10);

vp1 = nu(:, 1);
vp2 = nu(:, 2);
vp3 = nu(:, 3);
omega1 = rad2deg(nu(:, 4));
omega2 = rad2deg(nu(:, 5));
omega3 = rad2deg(nu(:, 6));
omegas1 = rad2deg(nu(:, 7));
omegas2 = rad2deg(nu(:, 8));
omegas3 = rad2deg(nu(:, 9));
xidot = nu(:, 10);

xi = y(:, 11);


% Compute Mechanical Energy
m = length(xi);
T = zeros(m, 1);
V = zeros(m, 1);

for i = 1 : m

    v = [vp1(i); vp2(i); vp3(i); omega1(i); omega2(i); omega3(i); omegas1(i); omegas2(i); omegas3(i); xidot(i)];
    Mt = buildMt(buildJP(xi(i)), buildSP(xi(i)));
    T(i) = 0.5 * v' * Mt * v;

    V(i) = 0.5 * Kd * xi(i)^2;

end

E = T + V;


%% Plot of the Results

figure('Name', "Body Angular Velocities")

hold on
plot(tspan, omega1)
plot(tspan, omega2)
plot(tspan, omega3)
title(strcat('Case', {' '}, choice, '- Plot of Body Angular Velocities'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\omega_i \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('\omega_1', '\omega_2', '\omega_3', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\omega-', choice, '.jpg'))
end


figure('Name', "Wheel Angular Velocities")

hold on
plot(tspan, omegas1)
plot(tspan, omegas2)
plot(tspan, omegas3, '--')
title(strcat('Case', {' '}, choice, '- Plot of Wheels Angular Velocities'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\omega s_i \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('\omega_{s1}', '\omega_{s2}', '\omega_{s3}', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\omegas-', choice, '.jpg'))
end


figure('Name', "Damper Mass Position")

hold on
plot(tspan, xi)
title(strcat('Case', {' '}, choice, '- Plot of Damper Mass Position'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\xi \ [m]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\xi-', choice, '.jpg'))
end


figure('Name', "Damper Mass Velocity")

hold on
plot(tspan, xidot)
title(strcat('Case', {' '}, choice, '- Plot of Damper Mass Velocity'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\dot{\xi} \ [m/s]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\xidot-', choice, '.jpg'))
end


figure('Name', "Mechanical Energy")

hold on
plot(tspan, E)
title(strcat('Case', {' '}, choice, '- Plot of Mechanical Energy'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$E \ [J]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\mechenergy-', choice, '.jpg'))
end
