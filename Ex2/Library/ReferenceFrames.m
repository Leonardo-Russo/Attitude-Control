%% Frames of References

close all
clear all
clc

B1 = eye(3);

theta1 = deg2rad(-30);

R21 = [1, 0, 0;...
       0, -cos(theta1), sin(theta1);...
       0, -sin(theta1), -cos(theta1)];

theta2 = deg2rad(30);

R31 = [1, 0, 0;...
       0, cos(theta2), -sin(theta2);...
       0, sin(theta2), cos(theta2)];

theta3 = deg2rad(-30);

R41 = [cos(theta3), 0, sin(theta3);...
       0, 1, 0;...
       -sin(theta3), 0, cos(theta3)];

theta4 = deg2rad(10);

R51 = [-cos(theta4), 0, -sin(theta4);...
       0, 1, 0;...
       sin(theta4), 0, -cos(theta4)];

B2 = R21' * B1;
B3 = R31' * B1;
B4 = R41' * B1;
B5 = R51' * B1;

indices = ['1', '2', '3'];

plot3(0,0,0,'g*')
title('Frames of Reference')
xlabel('x');
ylabel('y');
zlabel('z');
axis([-1, 1, -1, 1, -1, 1])
grid on

hold on

for i = 1 : 3
    quiver3(0, 0, 0, B1(i, 1), B1(i, 2), B1(i, 3), 'Color','b')
    text(B1(i, 1), B1(i, 2), B1(i, 3), strcat('B1,', indices(i)), 'Color','b')
end

for i = 1 : 3
    quiver3(0, 0, 0, B2(i, 1), B2(i, 2), B2(i, 3), 'Color', 'r')
    text(B2(i, 1), B2(i, 2), B2(i, 3), strcat('B2,', indices(i)), 'Color','r')
end

for i = 1 : 3
    quiver3(0, 0, 0, B3(i, 1), B3(i, 2), B3(i, 3), 'Color', 'r')
    text(B3(i, 1), B3(i, 2), B3(i, 3), strcat('B3,', indices(i)), 'Color','r')
end

for i = 1 : 3
    quiver3(0, 0, 0, B4(i, 1), B4(i, 2), B4(i, 3), 'Color', 'g')
    text(B4(i, 1), B4(i, 2), B4(i, 3), strcat('B4,', indices(i)), 'Color','g')
end

for i = 1 : 3
    quiver3(0, 0, 0, B5(i, 1), B5(i, 2), B5(i, 3), 'Color', 'g')
    text(B5(i, 1), B5(i, 2), B5(i, 3), strcat('B5,', indices(i)), 'Color','g')
end

% view(120, 0)