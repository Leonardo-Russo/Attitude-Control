%% Initialize

close all
clear all
clc

% addpath('Output\')

savefig = 1;

figure('Name','Reference Frames')

B1 = eye(3);
indices = ['1', '2', '3'];

plot3(0,0,0,'g*')
% title('Frames of Reference')
xlabel('x');
ylabel('y');
zlabel('z');
axis([-1, 1, -1, 1, -1, 1])
grid on

view(115, 30)

hold on
for i = 1 : 3
    quiver3(0, 0, 0, B1(i, 1), B1(i, 2), B1(i, 3), 'Color','b')
    text(B1(i, 1), B1(i, 2), B1(i, 3), strcat('$\hat{b_', indices(i), '}$'), 'Color','b', 'Interpreter','latex', 'FontSize', 15)
end

a4 = [sqrt(3)/3, sqrt(3)/3, sqrt(3)/3];
quiver3(0, 0, 0, a4(1), a4(2), a4(3), 'Color','b')
text(a4(1), a4(2), a4(3), strcat('$\hat{a_4}$'), 'Color','b', 'Interpreter','latex', 'FontSize', 15)

if savefig
    saveas(gcf, 'ReferenceFrames.jpg')
end