%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2
%% 
clear; clc;

%Initial conditions
omega_init = [0.05; 0.05; 0.05];
T = 2*pi/norm(omega_init);
%Integrate
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
tstart = 0; tint = 0.1; tend = 1*T/tint;
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

%% Ellipsoids
inertia_p = diag([85.075, 110.796, 120.2515]); %inertia in principal axes

%energy ellipsoid
T = 0.5* (omega_init(1)^2*inertia_p(1,1) + omega_init(2)^2*inertia_p(2,2) + omega_init(3)^2*inertia_p(3,3));
[x_enE, y_enE, z_enE] = ellipsoid(0,0,0, sqrt(2*T/inertia_p(1,1)), sqrt(2*T/inertia_p(2,2)), sqrt(2*T/inertia_p(3,3))); %energy ellipsoid

%momentum ellipsoid
L =  sqrt(omega_init(1)^2*inertia_p(1,1)^2 + omega_init(2)^2*inertia_p(2,2)^2 + omega_init(3)^2*inertia_p(3,3)^2);
[x_momE, y_momE, z_momE] = ellipsoid(0,0,0, sqrt(L^2/inertia_p(1,1)^2), sqrt(L^2/inertia_p(2,2)^2), sqrt(L^2/inertia_p(3,3)^2)); %energy ellipsoid


%% Validate Ellipsoid
%Energy
a = sqrt(2*T/inertia_p(1,1));
b = sqrt(2*T/inertia_p(2,2));
c = sqrt(2*T/inertia_p(3,3));
figure()
surface(x_enE, y_enE, z_enE, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);
hold on;
plot3(linspace(0,a), 0, 0,'.', 'Color','black', 'LineWidth', 2)
plot3(0,linspace(0,b), 0, '.', 'Color','black', 'LineWidth', 2)
plot3(0, 0, linspace(0,c),'.', 'Color','black', 'LineWidth', 2)
view(3)
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
zlabel('\omega_z [1/s^2]')
grid on;
axis equal;
title('Energy Ellipsoid with semi-major axes')
hold off;

%Momentum
a = sqrt(L^2/inertia_p(1,1)^2);
b = sqrt(L^2/inertia_p(2,2)^2);
c = sqrt(L^2/inertia_p(3,3)^2);
figure()
surface(x_momE, y_momE, z_momE, 'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);
hold on;
plot3(linspace(0,a), 0, 0,'.', 'Color','black', 'LineWidth', 2)
plot3(0,linspace(0,b), 0, '.', 'Color','black', 'LineWidth', 2)
plot3(0, 0, linspace(0,c),'.', 'Color','black', 'LineWidth', 2)
view(3)
grid on;
axis equal;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
zlabel('\omega_z [1/s^2]')
title('Momentum Ellipsoid with semi-major axes')

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

% Side views
figure();
subplot(2,2,1);
plot(omega_out(:,1), omega_out(:,2))
axis equal;
grid on;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
legend('Polhode in XY-plane')

subplot(2,2,2);
plot(omega_out(:,1), omega_out(:,3))
axis equal;
grid on;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_z [1/s^2]')
legend('Polhode in XZ-plane')

subplot(2,2,3);
plot(omega_out(:,2), omega_out(:,3))
axis equal;
grid on;
xlabel('\omega_y [1/s^2]')
ylabel('\omega_z [1/s^2]')
legend('Polhode in YZ-plane')
%% return omega_dot (zero torque)
function [omegadot] = omega_dot(t, omega)
    inertia_p = diag([85.075, 110.796, 120.2515]);
    M = zeros(3,1); %no torques
    %Euler eqns in principal axes
    omegadot = [ 1/inertia_p(1,1) * (M(1) - (inertia_p(3,3) - inertia_p(2,2)) * omega(2)*omega(3));
                 1/inertia_p(2,2) * (M(2) - (inertia_p(1,1) - inertia_p(3,3)) * omega(3)*omega(1));
                 1/inertia_p(3,3) * (M(3) - (inertia_p(2,2) - inertia_p(1,1)) * omega(1)*omega(2))];
end