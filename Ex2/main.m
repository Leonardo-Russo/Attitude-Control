%% Homework 2 - Leonardo Russo 2015563

close all
clear all
clc

addpath('Library\')

%% Introduction to Known Quantities

R0 = 7000;      % km
mu = 398600.4415;       % km^3/s^2

global mB1 mBj R L r lA lB K1 K2

mB1 = 2000;     % kg
mBj = 100;      % kg

R = 1;      % m
L = 4;      % m
r = 0.05;   % m

K1 = 100;   % Nm
K2 = 50;    % Nms


%% Initial Conditions

% Define the Reference Frames
theta1_0 = deg2rad(-30);    % rad
theta2_0 = deg2rad(30);     % rad
theta3_0 = deg2rad(-30);    % rad
theta4_0 = deg2rad(10);     % rad
thetas_0 = [theta1_0 theta2_0 theta3_0 theta4_0]';

vp1_0 = [0 sqrt(mu/R0) 0]' * 1e3;   % m/s
rp1_0 = [R0 0 0]' * 1e3;            % m
omega1_0 = [1 0.1 0.1]';        % rad/s

sigma1_0 = -0.2;    % rad/s
sigma2_0 = 0.2;     % rad/s
sigma3_0 = 0.1;     % rad/s
sigma4_0 = -0.3;    % rad/s
sigmas_0 = [sigma1_0 sigma2_0 sigma3_0 sigma4_0]';

R1N_0 = eye(3);     % Initial Attitude of B1 wrt N
[q10, q1] = C2q(R1N_0);     % convert Attitude in Quaternions

%% Options

t0 = 0;     % sec
lA = 2;     % m

% choice = '2';
% savechoice = '0';

choice = char(string(input('Please, select the case:\n1. Case 1\n2. Case 2\n')));
clc

savechoice = char(string(input('Do you wish to save a copy of the graphs locally?\n1. Yes\n2. No\n')));
clc

switch choice

    case '1'
        tf = 1800;  % sec
        lB = 9;     % m

    case '2'
        tf = 7200;  % sec
        lB = 3;     % m
end


%% Numerical Integration

% Define the options for the Integration
Tol0 = 1e-11;
Tol1 = 1e-13;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

% Define Initial State
x0 = [q10; q1];
X0 = [omega1_0; sigmas_0; vp1_0; rp1_0; thetas_0; x0];

% Perform the Integration
[tspan, X] = ode113('DynamicalModel', [t0, tf], X0, options);

% Extract the Results
omega1 = rad2deg(X(:, 1:3));
sigmas = rad2deg(X(:, 4:7));
vp1 = X(:, 8:10) * 1e-3;
thetas = rad2deg(X(:, 14:17));


%% Plot of the Results

figure('Name', "Central Body Angular Velocity")

hold on
plot(tspan, omega1(:, 1))
plot(tspan, omega1(:, 2))
plot(tspan, omega1(:, 3))
% title(strcat('Case', {' '}, choice, '- Plot of Central Body Angular Velocity'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\omega_1} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('\omega_{1,1}', '\omega_{1,2}', '\omega_{1,3}', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\omega1-', choice, '.jpg'))
end


figure('Name', "Central Body Linear Velocity")

hold on
plot(tspan, vp1(:, 1))
plot(tspan, vp1(:, 2))
plot(tspan, vp1(:, 3))
% title(strcat('Case', {' '}, choice, '- Plot of Central Body Velocity'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{v_{p1}} \ [km/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('v_{p1,1}', 'v_{p1,2}', 'v_{p1,3}', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\vp1-', choice, '.jpg'))
end


figure('Name', "Joint Angular Rate 1")

hold on
plot(tspan, sigmas(:, 1))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\sigma_{1} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\sigma1-', choice, '.jpg'))
end


figure('Name', "Joint Angular Rate 2")

hold on
plot(tspan, sigmas(:, 2))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\sigma_{2} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\sigma2-', choice, '.jpg'))
end


figure('Name', "Joint Angular Rate 3")

hold on
plot(tspan, sigmas(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\sigma_{3} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\sigma3-', choice, '.jpg'))
end


figure('Name', "Joint Angular Rate 4")

hold on
plot(tspan, sigmas(:, 4))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\sigma_{4} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\sigma4-', choice, '.jpg'))
end


figure('Name', "Joint Angular Rates 1-2")

hold on
plot(tspan, sigmas(:, 1))
plot(tspan, sigmas(:, 2))
% title(strcat('Case', {' '}, choice, '- Plot of \sigma_1 and \sigma_2'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\sigma_{i}} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('\sigma_1', '\sigma_2', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\sigmas12-', choice, '.jpg'))
end


figure('Name', "Joint Angular Rates 3-4")

hold on
plot(tspan, sigmas(:, 3))
plot(tspan, sigmas(:, 4))
% title(strcat('Case', {' '}, choice, '- Plot of \sigma_3 and \sigma_4'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\sigma_{i}} \ [deg/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('\sigma_3', '\sigma_4', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\sigmas34-', choice, '.jpg'))
end


figure('Name', "Joint Angle 1")

hold on
plot(tspan, thetas(:, 1))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\theta_{1} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\theta1-', choice, '.jpg'))
end


figure('Name', "Joint Angle 2")

hold on
plot(tspan, thetas(:, 2))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\theta_{2} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\theta2-', choice, '.jpg'))
end


figure('Name', "Joint Angle 3")

hold on
plot(tspan, thetas(:, 3))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\theta_{3} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\theta3-', choice, '.jpg'))
end


figure('Name', "Joint Angle 4")

hold on
plot(tspan, thetas(:, 4))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\theta_{4} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\theta4-', choice, '.jpg'))
end


figure('Name', "Joint Angles 1-2")

hold on
plot(tspan, thetas(:, 1))
plot(tspan, thetas(:, 2))
% title(strcat('Case', {' '}, choice, '- Plot of \theta_1 and \theta_2'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\theta_{i}} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
legend('\theta_1', '\theta_2', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\thetas12-', choice, '.jpg'))
end


figure('Name', "Joint Angles 3-4")

hold on
plot(tspan, thetas(:, 3))
plot(tspan, thetas(:, 4))
% title(strcat('Case', {' '}, choice, '- Plot of \theta_3 and \theta_4'))
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\underline{\theta_{i}} \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
legend('\theta_3', '\theta_4', 'Location','best')
if savechoice == '1'
    saveas(gcf, strcat('Output\thetas34-', choice, '.jpg'))
end

