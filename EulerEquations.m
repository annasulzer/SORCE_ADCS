%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2
%% 
clear; clc; close all;

%Initial conditions
 omega_init = [0.1; 0.15; 0.08];
% omega_init = [0.25; 0.001; 0.001];
% omega_init = [0.001; 0.25; 0.001];
% omega_init = [0.001; 0.001; 0.25];

time = 2*pi/norm(omega_init);
%Integrate
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
tstart = 0; tint = 0.1; tend = 10*time/tint;
[t_out, omega_out] = ode113(@omega_dot,[tstart:tint:tend]', omega_init, options);
t_out = t_out * tint; %to make seconds


% Plot
figure()
subplot(3,1,1)
plot(t_out, omega_out(:, 1))
xlabel('t [s]')
ylabel('\omega_x [rad/s^2]')
title('Angular Velocity over time')
subplot(3,1,2)
plot(t_out, omega_out(:, 2))
xlabel('t [s]')
ylabel('\omega_y [rad/s^2]')
subplot(3,1,3)
plot(t_out, omega_out(:, 3))
xlabel('t [s]')
ylabel('\omega_z [rad/s^2]')

saveas(gcf, 'ang_vel.pdf');
%% Ellipsoids
[R_princ,inertia_p] = inertia(); %inertia in principal axes

%energy ellipsoid
T = 0.5* (omega_init(1)^2*inertia_p(1,1) + omega_init(2)^2*inertia_p(2,2) + omega_init(3)^2*inertia_p(3,3));
[x_enE, y_enE, z_enE] = ellipsoid(0,0,0, sqrt(2*T/inertia_p(1,1)), sqrt(2*T/inertia_p(2,2)), sqrt(2*T/inertia_p(3,3))); %energy ellipsoid

%Semi Major Axes
aE = sqrt(2*T/inertia_p(1,1));
bE = sqrt(2*T/inertia_p(2,2));
cE = sqrt(2*T/inertia_p(3,3));

%momentum ellipsoid
L =  sqrt(omega_init(1)^2*inertia_p(1,1)^2 + omega_init(2)^2*inertia_p(2,2)^2 + omega_init(3)^2*inertia_p(3,3)^2);
[x_momE, y_momE, z_momE] = ellipsoid(0,0,0, sqrt(L^2/inertia_p(1,1)^2), sqrt(L^2/inertia_p(2,2)^2), sqrt(L^2/inertia_p(3,3)^2)); %energy ellipsoid

%Momentum
aM = sqrt(L^2/inertia_p(1,1)^2);
bM = sqrt(L^2/inertia_p(2,2)^2);
cM = sqrt(L^2/inertia_p(3,3)^2);

%Check that it's real
disp(L^2/(2*T)); %to check if real


% Analytical Solutions
%yz plane
%Reference ellipse
a_yz = sqrt((L^2-2*T*inertia_p(1,1))/((inertia_p(2,2) - inertia_p(1,1))*inertia_p(2,2)));
b_yz = sqrt((L^2-2*T*inertia_p(1,1))/((inertia_p(3,3) - inertia_p(1,1))*inertia_p(3,3)));

%xz plane
%Reference hyperbola
a_xz = sqrt((L^2 -2*T*inertia_p(2,2))/(((inertia_p(1,1) - inertia_p(2,2)) * inertia_p(1,1))));
b_xz = sqrt((L^2 -2*T*inertia_p(2,2))/(((inertia_p(3,3) - inertia_p(2,2)) * inertia_p(3,3))));

%xy plane
%Reference ellipse
a_xy = sqrt((L^2 -2*T*inertia_p(3,3))/(((inertia_p(1,1) - inertia_p(3,3)) * inertia_p(1,1))));
b_xy = sqrt((L^2 -2*T*inertia_p(3,3))/(((inertia_p(2,2) - inertia_p(3,3)) * inertia_p(2,2))));



%% Plot & Validate Ellipsoid
%Energy
figure()
surface(x_enE, y_enE, z_enE, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);
hold on;
plot3(linspace(0,aE), 0, 0,'.', 'Color','black', 'LineWidth', 2)
plot3(0,linspace(0,bE), 0, '.', 'Color','black', 'LineWidth', 2)
plot3(0, 0, linspace(0,cE),'.', 'Color','black', 'LineWidth', 2)
view(3)
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
zlabel('\omega_z [1/s^2]')
grid on;
axis equal;
title('Energy Ellipsoid with semi-major axes')
hold off;
saveas(gcf, 'EnEllipsoid.pdf');

%Momentum
figure()
surface(x_momE, y_momE, z_momE, 'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);
hold on;
plot3(linspace(0,aM), 0, 0,'.', 'Color','black', 'LineWidth', 2)
plot3(0,linspace(0,bM), 0, '.', 'Color','black', 'LineWidth', 2)
plot3(0, 0, linspace(0,cM),'.', 'Color','black', 'LineWidth', 2)
view(3)
grid on;
axis equal;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
zlabel('\omega_z [1/s^2]')
title('Momentum Ellipsoid with semi-major axes')
%saveas(gcf, 'MomEllipsoid.pdf');
%% Plot
%Ellipsoid
figure()
surface(x_enE, y_enE, z_enE, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);
hold on;
surface(x_momE, y_momE, z_momE, 'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);

%Polhode
plot3(omega_out(:,1), omega_out(:,2), omega_out(:, 3), 'Color', 'black', 'LineWidth', 2)

axis equal;
grid on;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
zlabel('\omega_z [1/s^2]')
legend('Energy Ellipsoid', 'Momentum Ellipsoid', 'Polhode')
view(3)
%saveas(gcf, 'Polhode.pdf');
%% Plot Side views
%yz plane
figure();
subplot(2,2,1);
hold on;
plot(omega_out(:,2), omega_out(:,3),'Color', 'red', 'LineWidth',2)% data
fplot(@(y) sqrt(b_yz^2 * (1 - y^2 / a_yz^2)), '--','LineWidth',2, 'Color','blue'); % Upper half of ellipse
fplot(@(y) -sqrt(b_yz^2 * (1 - y^2 / a_yz^2)), '--', 'LineWidth',2, 'Color','blue'); % Upper half of ellipse
axis equal;
grid on;
xlabel('\omega_y [1/s^2]')
ylabel('\omega_z [1/s^2]')
title('Polhode in YZ-plane')

%xz plane
subplot(2,2,2);
hold on;
plot(omega_out(:,1), omega_out(:,3),'Color', 'red', 'LineWidth',2)% data
fplot(@(x) sqrt(b_xz^2 * (1 - x^2 / a_xz^2)),'--', 'LineWidth',2, 'Color','blue'); % half of hyperbola
fplot(@(x) -sqrt(b_xz^2 * (1 - x^2 / a_xz^2)), '--','LineWidth',2, 'Color','blue'); % half of hyperbola
axis equal;
grid on;
xlim([-1.5*norm(max(omega_out)), 1.5*norm(max(omega_out))]);
ylim([-1.5*norm(max(omega_out)), 1.5*norm(max(omega_out))]);
xlabel('\omega_x [1/s^2]')
ylabel('\omega_z [1/s^2]')
title('Polhode in XZ-plane')
hold off;

%xy plane
subplot(2,2,3);
hold on;
plot(omega_out(:,1), omega_out(:,2),'Color', 'red', 'LineWidth',2)% data)
fplot(@(x) sqrt(b_xy^2 * (1 - x^2 / a_xy^2)), '--','LineWidth',2, 'Color','blue'); % Upper half of ellipse
fplot(@(x) -sqrt(b_xy^2 * (1 - x^2 / a_xy^2)),'--', 'LineWidth',2, 'Color','blue'); % Upper half of ellipse
axis equal;
grid on;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
title('Polhode in XY-plane')
hold off;
legend('Simulation Data', 'Analytical Solution', '')
%saveas(gcf, 'Polhode2D.pdf');

%% return omega_dot (zero torque)
function [omegadot] = omega_dot(t, omega)
    [R_princ,inertia_p] = inertia();
    M = zeros(3,1); %no torques
    %Euler eqns in principal axes
    omegadot = [ 1/inertia_p(1,1) * (M(1) - (inertia_p(3,3) - inertia_p(2,2)) * omega(2)*omega(3));
                 1/inertia_p(2,2) * (M(2) - (inertia_p(1,1) - inertia_p(3,3)) * omega(3)*omega(1));
                 1/inertia_p(3,3) * (M(3) - (inertia_p(2,2) - inertia_p(1,1)) * omega(1)*omega(2))];
end