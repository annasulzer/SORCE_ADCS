%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2 and PSET3
%  For Polhode and Ellipsoid Plotting run EllipsoidAndPolhode.m
%% 
clear; clc; 
close all;
%%
%Initial conditions
[state_ECI_init, T_orbit] = OrbitPropagation();
omega_init = [0.0; 0.0; 0.0];

%IC based on EULER angle
att_init = [0, 90, 0]; %313
Rot = rotz(-att_init(3))*rotx(-att_init(2))*rotz(-att_init(1));

%IC aligned with RTN frame
% [R, T, N] = RTN_frame_inertial(state_ECI_init');
% Rot = [R',T', N']';

DCM_initial = Rot;

[R_princ,inertia_p] = inertia(); %inertia in principal axes
%inertia_p = [85.075, 0, 0; 0, 85.075, 0; 0, 0, 120.2515]; %axial symmetric

%Integration settings
absTol= 1e-10;
relTol = 1e-6;
%time = 2*pi/norm(omega_init);
tstart = 0; tend = T_orbit;

%% Simulate
out = sim("main");
%% Read out simulation data
%attitude propagation
t_out = out.tout;
omega_out = out.omega.Data(:,:)';
L_out = (inertia_p * omega_out')';

%attitude representation SWITCH
quat_out = out.quaternions.Data(:, :)';
euler_out = out.euler.Data(:, :)';

omega_out_inertial =  out.omega_inertial.Data(:, :)';
L_out_inertial =  out.L_inertial.Data(:, :)';
DCM_out = out.DCM.Data(:,:,:);

%orbit propagation
state_out = out.orbit_state.Data(:,:)';
torque_out = out.M.Data(:,:)';

%% Get different coordinate frames with respect to inertial frame
Rot = DCM_out;
[Xp, Yp, Zp, Xb, Yb, Zb] = principal_body_frame_inertial(Rot, R_princ);
[R, T, N] = RTN_frame_inertial(state_out);

%% Plotting PS2
% Plot omega over time
% figure()
% subplot(3,1,1)
% plot(t_out, omega_out(:, 1))
% xlabel('t [s]')
% ylabel('\omega_x [rad/s]')
% title('Angular Velocity over time from the numerical integration')
% subplot(3,1,2)
% plot(t_out, omega_out(:, 2))
% xlabel('t [s]')
% ylabel('\omega_y [rad/s]')
% subplot(3,1,3)
% plot(t_out, omega_out(:, 3))
% xlabel('t [s]')
% ylabel('\omega_z [rad/s]')
% 
% % Analytical Solution for axially symmetric satellite
% %Angular Velocity
% omega_analyt=zeros(length(t_out), 3);
% lambda = omega_init(3)*(inertia_p(3,3) - inertia_p(1,1))/inertia_p(1,1);
% theta_zero = atan2(omega_init(2), omega_init(1));
% for i = 1:length(t_out)
%     omega_analyt(i, 3) = omega_init(3); %omega_z 
%     omega_analyt(i, 1) = norm(omega_init(1:2)) * cos(lambda*t_out(i) + theta_zero); %omega_x = om_xy * cos(lambda*t)
%     omega_analyt(i, 2) = norm(omega_init(1:2)) * sin(lambda*t_out(i) + theta_zero); %omega_x = om_xy * sin(lambda*t)
% end
% %Errors btw analytical and numerical
% errors = omega_out - omega_analyt;
% 
% % Plot over time
% figure()
% subplot(3,1,1)
% plot(t_out, omega_analyt(:, 1))
% xlabel('t [s]')
% ylabel('\omega_x [rad/s]')
% title('Angular Velocity over time from the analytical solution')
% subplot(3,1,2)
% plot(t_out, omega_analyt(:, 2))
% xlabel('t [s]')
% ylabel('\omega_y [rad/s]')
% subplot(3,1,3)
% plot(t_out, omega_analyt(:, 3))
% xlabel('t [s]')
% ylabel('\omega_z [rad/s]')
% 
% %Plot Error over time
% figure()
% subplot(3,1,1)
% hold on;
% plot(t_out, errors(:, 1), 'red')
% xlabel('t [s]')
% ylabel('\Delta\omega_x [rad/s]')
% title('Errors between numerical solution and analytical solution over time')
% subplot(3,1,2)
% hold on;
% plot(t_out, errors(:, 2), 'red')
% xlabel('t [s]')
% ylabel('\Delta\omega_y [rad/s]')
% subplot(3,1,3)
% hold on;
% plot(t_out, errors(:, 3), 'red')
% xlabel('t [s]')
% ylabel('\Delta\omega_z [rad/s]')
% 
% %% Plotting Angular Momentum Vector
% % Analytical Solution for axially symmetric satellite
% %Angular Momentum
% L_analyt=zeros(length(t_out), 3);
% L_analyt_init = inertia_p * omega_init;
% for i = 1:length(t_out)
%     L_analyt(i, 3) = L_analyt_init(3); %omega_z 
%     L_analyt(i, 1) = norm(L_analyt_init(1:2)) * cos(lambda*t_out(i) + theta_zero); %omega_x = om_xy * cos(lambda*t)
%     L_analyt(i, 2) = norm(L_analyt_init(1:2)) * sin(lambda*t_out(i) + theta_zero); %omega_x = om_xy * sin(lambda*t)
% end
% 
% 
% %Plot Angular Momentum
% figure()
% subplot(3,1,1)
% hold on;
% plot(t_out, L_out(:, 1), 'black')
% plot(t_out, L_analyt(:, 1),'--', 'Color','green')
% xlabel('t [s]')
% ylabel('L_x [kg*m^2/s]')
% title('Angular momentum vector of numerical solution and analytical solution over time')
% subplot(3,1,2)
% hold on;
% plot(t_out, L_out(:, 2), 'black')
% plot(t_out, L_analyt(:, 2), '--','Color','green')
% xlabel('t [s]')
% ylabel('L_y [kg*m^2/s]')
% subplot(3,1,3)
% hold on;
% 
% plot(t_out, L_out(:, 3), 'black')
% plot(t_out, L_analyt(:, 3),'--', 'Color','green')
% xlabel('t [s]')
% ylabel('L_z [kg*m^2/s]')

%% Plot Attitude Representations PS3
% Quaternions over time
% figure()
% subplot(4,1,1)
% hold on;
% plot(t_out, quat_out(:, 4), 'blue')
% xlabel('t [s]')
% ylabel('q4')
% title('Quaternions over time')
% subplot(4,1,2)
% hold on;
% plot(t_out, quat_out(:, 1), 'blue')
% xlabel('t [s]')
% ylabel('q1')
% subplot(4,1,3)
% hold on;
% plot(t_out, quat_out(:, 2), 'blue')
% xlabel('t [s]')
% ylabel('q2')
% subplot(4,1,4)
% hold on;
% plot(t_out, quat_out(:, 3), 'blue')
% xlabel('t [s]')
% ylabel('q3')
% 
% % Euler Angles over time
% figure()
% subplot(3,1,1)
% hold on;
% plot(t_out, rad2deg(euler_out(:, 1)), 'blue')
% xlabel('t [s]')
% ylabel('\phi [deg]')
% title('Euler angles over time')
% subplot(3,1,2)
% hold on;
% plot(t_out, rad2deg(euler_out(:, 2)), 'blue')
% xlabel('t [s]')
% ylabel('\theta [deg]')
% subplot(3,1,3)
% hold on;
% plot(t_out, rad2deg(euler_out(:, 3)), 'blue')
% xlabel('t [s]')
% ylabel('\psi [deg]')
% 
% %% Check Attitude Representation
% % Angular Momentum Vector
% figure()
% subplot(3,1,1)
% hold on;
% plot(t_out, L_inertial_euler_out(:, 1), 'blue')
% plot(t_out, L_inertial_quat_out(:, 1), '--', 'Color','red')
% xlabel('t [s]')
% ylabel('L_x [kg*m^2/s]')
% title('Inertial angular momentum vector over time from the euler angles (blue) and quaternions (red)')
% subplot(3,1,2)
% hold on;
% plot(t_out, L_inertial_euler_out(:, 2), 'blue')
% plot(t_out, L_inertial_quat_out(:, 2), '--', 'Color','red')
% xlabel('t [s]')
% ylabel('L_y [kg*m^2/s]')
% subplot(3,1,3)
% hold on;
% plot(t_out, L_inertial_euler_out(:, 3), 'blue')
% plot(t_out, L_inertial_quat_out(:, 3), '--', 'Color','red')
% xlabel('t [s]')
% ylabel('L_z [kg*m^2/s]')
% legend('Euler Angles', 'Quaternions')
% 
% %%
% % Plot Omega over time
% figure()
% subplot(3,1,1)
% hold on;
% plot(t_out, omega_inertial_euler_out(:, 1), 'blue')
% plot(t_out, omega_inertial_quat_out(:, 1), '--', 'Color','red')
% xlabel('t [s]')
% ylabel('\omega_x [rad/s]')
% title('Inertial angular velocity over time from the euler angles (blue) and quaternions (red)')
% subplot(3,1,2)
% hold on;
% plot(t_out, omega_inertial_euler_out(:, 2), 'blue')
% plot(t_out, omega_inertial_quat_out(:, 2), '--', 'Color','red')
% xlabel('t [s]')
% ylabel('\omega_y [rad/s]')
% subplot(3,1,3)
% hold on;
% plot(t_out, omega_inertial_euler_out(:, 3), 'blue')
% plot(t_out, omega_inertial_quat_out(:, 3), '--', 'Color','red')
% xlabel('t [s]')
% ylabel('\omega_z [rad/s]')
% legend('Euler Angles', 'Quaternions')
% 
% %Herpolhode Quaternions
% figure()
% hold on;
% plot3(omega_inertial_quat_out(:,1), omega_inertial_quat_out(:,2), omega_inertial_quat_out(:, 3), 'Color', 'black', 'LineWidth', 1)
% quiver3(0,0,0, L_inertial_quat_out(1,1)/50, L_inertial_quat_out(1,2)/50, L_inertial_quat_out(1,3)/50,'LineWidth', 1)
% axis equal;
% grid on;
% xlabel('x')
% ylabel('y')
% zlabel('z')
% legend('Herpolhode', 'Angular Momentum Vector (scaled)')
% title('Herpolhode and Angular Momentum Vector based on Quaternions')
% view(3)
% 
% %Herpolhode Euler
% figure()
% hold on;
% plot3(omega_inertial_euler_out(:,1), omega_inertial_euler_out(:,2), omega_inertial_euler_out(:, 3), 'Color', 'black', 'LineWidth', 1)
% quiver3(0,0,0, L_inertial_euler_out(1,1)/50, L_inertial_euler_out(1,2)/50, L_inertial_euler_out(1,3)/50,'LineWidth', 1)
% axis equal;
% grid on;
% xlabel('x')
% ylabel('y')
% zlabel('z')
% legend('Herpolhode', 'Angular Momentum Vector (scaled)')
% title('Herpolhode and Angular Momentum Vector based on Euler Angles')
% view(3)
% 
%% Plotting PSET4
%Euler Angles over time
figure()
hold on;
plot(t_out, rad2deg(unwrap(euler_out(:, 1))), LineWidth=2)
plot(t_out, rad2deg(unwrap(euler_out(:, 2))), LineWidth=2)
plot(t_out, rad2deg(unwrap(euler_out(:, 3))), LineWidth=2)
xlabel('t [s]')
ylabel('Angle [deg]')
legend('\phi', '\theta', '\psi')
title('Euler angles over time (213 sequence)')

%Plot Omega over time
figure()
hold on;
plot(t_out, omega_out_inertial(:, 1), LineWidth=2)
plot(t_out, omega_out_inertial(:, 2), LineWidth=2)
plot(t_out, omega_out_inertial(:, 3), LineWidth=2)
xlabel('t [s]')
ylabel('\omega [rad/s]')
legend('X', 'Y', 'Z')
title('Inertial angular velocity over time')


%% Problem 4
% Verify magnitude
%Plot Torque over time
figure()
hold on;
plot(t_out, torque_out(:, 1), LineWidth=2)
plot(t_out, torque_out(:, 2), '--',LineWidth=2)
plot(t_out, torque_out(:, 3),LineWidth=2)
plot(t_out, sqrt(torque_out(:, 1).^2 + torque_out(:, 2).^2 + torque_out(:, 3).^2),'--', LineWidth=2)
xlabel('t [s]')
ylabel('M [Nm]')
legend('M_x', 'M_y', 'M_z', 'M_{tot}')
title('Torque over time (one orbit)')

