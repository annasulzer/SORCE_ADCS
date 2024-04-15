%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2
%% 
clear; clc;

%Initial conditions
omega_init = [0; 1; 0.001];
rps = norm(omega_init)/2*pi; %rotations per second

%Integrate
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
tstart = 0; tint = 0.01; tend = 1000/rps;
[t_out, omega_out] = ode113(@omega_dot,[tstart:tint:tend]', omega_init, options);

%% Ellipsoids
inertia_p = diag([85.075, 110.796, 120.2515]); %inertia in principal axes

%energy ellipsoid
T = 0.5* (omega_init(1)^2*inertia_p(1,1) + omega_init(2)^2*inertia_p(2,2) + omega_init(3)^2*inertia_p(3,3));
[x_enE, y_enE, z_enE] = ellipsoid(0,0,0, sqrt(2*T/inertia_p(1,1)), sqrt(2*T/inertia_p(2,2)), sqrt(2*T/inertia_p(3,3))); %energy ellipsoid

%momentum ellipsoid
L =  sqrt(omega_init(1)^2*inertia_p(1,1)^2 + omega_init(2)^2*inertia_p(2,2)^2 + omega_init(3)^2*inertia_p(3,3)^2);
[x_momE, y_momE, z_momE] = ellipsoid(0,0,0, sqrt(L^2/inertia_p(1,1)^2), sqrt(L^2/inertia_p(2,2)^2), sqrt(L^2/inertia_p(3,3)^2)); %energy ellipsoid


% %%
% figure()
% hold on;
% fsurf(@(omx,omy)  sqrt(2*T/inertia_p(3,3)*(1 - omx^2/(2*T/inertia_p(1,1)) - omy^2/(2*T/inertia_p(2,2)))), 'FaceColor', 'blue', 'EdgeColor', 'blue', 'FaceAlpha',.2)
% fsurf(@(omx,omy) -sqrt(2*T/inertia_p(3,3)*(1 - omx^2/(2*T/inertia_p(1,1)) - omy^2/(2*T/inertia_p(2,2)))), 'FaceColor', 'blue', 'EdgeColor', 'blue', 'FaceAlpha',.2)
% axis equal;
% 
% fsurf(@(omx,omy)  sqrt(L^2/inertia_p(3,3)^2*(1 - omx^2/(L^2/inertia_p(1,1)^2) - omy^2/(L^2/inertia_p(2,2)^2))), [-norm(omega_init) norm(omega_init) -norm(omega_init) norm(omega_init)], 'FaceColor', 'red', 'EdgeColor', 'red', 'FaceAlpha',.2)
% fsurf(@(omx,omy) -sqrt(L^2/inertia_p(3,3)^2*(1 - omx^2/(L^2/inertia_p(1,1)^2) - omy^2/(L^2/inertia_p(2,2)^2))), [-norm(omega_init) norm(omega_init) -norm(omega_init) norm(omega_init)], 'FaceColor', 'red', 'EdgeColor', 'red', 'FaceAlpha',.2)
% xlabel('\omega_x [1/s^2]')
% ylabel('\omega_y [1/s^2]')
% zlabel('\omega_z [1/s^2]')
% legend('Energy Ellipsoid','', 'Momentum Ellipsoid')
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