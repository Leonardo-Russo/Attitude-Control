function DrawTraj_3D(r_vect, xe)
% Description:
% Create a 3D Plot of the propagated orbit and save it into
% a pdf file.

X = r_vect(:, 1);
Y = r_vect(:, 2);
Z = r_vect(:, 3);

qe0 = xe(:, 1);
qe = xe(:, 2:4);
qe_tol = 0.005;

M = length(X);

X_ss = [];
X_tr = [];
Y_ss = [];
Y_tr = [];
Z_ss = [];
Z_tr = [];

for j = 1 : M
    
    if (qe(j, 1) > qe_tol || qe(j, 2) > qe_tol || qe(j, 3) > qe_tol) && j > floor(M/3)
        X_ss = [X_ss; X(j)];
        Y_ss = [Y_ss; Y(j)];
        Z_ss = [Z_ss; Z(j)];
    else
        X_tr = [X_tr; X(j)];
        Y_tr = [Y_tr; Y(j)];
        Z_tr = [Z_tr; Z(j)];
    end

end

fig_title = strcat("Orbit of the SpaceCraft");

figure('Name', fig_title, 'NumberTitle', 'off')

[x,y,z]=sphere;

I = imread('terra.jpg');

surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct')

hold on
plot3(X_ss,Y_ss,Z_ss,'Color','#ff7403', 'Marker', '.','MarkerSize', 10, 'Linestyle', 'none')
plot3(X_tr,Y_tr,Z_tr,'Color','#34c1e0', 'Marker', '.','MarkerSize', 5, 'Linestyle', 'none')
plot3(X(1), Y(1), Z(1), 'Color', '#eb2a71', 'Marker', '.', 'MarkerSize',15)
text(X(1)+500, Y(1)+500, Z(1)-500, '$\underline{x_0}$', 'Interpreter','latex', 'FontSize', 15, 'Color','#eb2a71')
plot3(0,0,0,'g*')
hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([130,25])


end