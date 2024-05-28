%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2 and PSET3
%  For Polhode and Ellipsoid Plotting run EllipsoidAndPolhode.m
%% 
clear; clc; 
close all;
%%
[R_princ,inertia_p] = inertia(); %inertia in principal axes

%Initial conditions
[state_ECI_init, T_orbit, n] = OrbitPropagation();
[t0_MJD_sun, a_sun, Om_sun, e_sun, om_sun, i_sun, M0_sun, n_sun] = OrbitPropagation_Sun();
%omega_init = R_princ * [0.0001; 0.0001; 0.1];
%omega_init =[n*0.1; n*0.1; n];

%inertia_p = [85.075, 0, 0; 0, 85.075, 0; 0, 0, 120.2515]; %axial symmetric

%Initial conditions
[state_ECI_init, T_orbit, n] = OrbitPropagation();
% omega_init = R_princ * [n; 0; 0];
omega_init = [0; 0; 0.0];

% Sun IC
UT1 = [1,25,2004,00];
JD_init = 2453029.5;
MJD_init = JD_init - 2400000.5;

D_init = JD_init - 2451545.0;
theta = UT1_to_theta(UT1);

%IC based on EULER angle
att_init = [0, 0, 0]; %313
% Rot2 = roty(-att_init(1))*rotx(-att_init(2))*rotz(-att_init(3));
% 
% %IC aligned with RTN frame
% [R, T, N] = RTN_frame_inertial(state_ECI_init');
% Rot = [R',T', N']';
% 
% DCM_initial = Rot;
% DCM_initial = R_princ' * Rot * Rot2;

DCM_initial = targetDCM([26321453.5527815,	-132781955.130633,	-57571626.5531097]', R_princ); %rinitial state sun

%Integration settings
eps = 1e-10;
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
Meas_sun = out.simout.Data(1:3,:)';
Meas_sunV = out.simout1.Data(1:3,:)';

%attitude representation SWITCH
quat_out = out.quaternions.Data(:, :)';
euler_out = out.euler.Data(:, :)';

omega_out_inertial =  out.omega_inertial.Data(:, :)';
L_out_inertial =  out.L_inertial.Data(:, :)';
DCM_out = out.DCM.Data;

%orbit propagation
state_out = out.orbit_state.Data(:,:)';
state_sun_out = out.state_sun_ECI.Data(:,:)';
M_grav_out = out.M_grav.Data(:,:)';
M_mag_out = out.M_mag.Data(:,:)';
M_SRP_out = out.M_SRP.Data(:,:)';
M_aero_out = out.M_aero.Data(:,:)';
M_ALL_out = out.M_ALL.Data(:,:)';

%State Estimation
DCM_estimated_det_out = out.DCM_estimated_det.Data;
quat_estimated_Q_out = out.quat_estimated_Q.Data(:, :)';
quat_estimated_kin_out = out.quat_estimated_kin.Data(:, :)';


eclipse_condition = out.eclipse.Data();
%% check that norm of one is still given
figure()
x = zeros(length(Meas_sun), 3);
y = zeros(length(Meas_sun), 3);
hold on;
for i = 1:length(Meas_sun)
    x(i, :) = (Meas_sun(i, :));
    y(i, :) = DCM_out(:,:,i)*Meas_sunV(i, :)';
end
plot(t_out, x(:, 1))
plot(t_out, y(:, 1), LineWidth=2)
plot(t_out, x(:, 2))
plot(t_out, y(:, 2), LineWidth=2)
plot(t_out, x(:, 3))
plot(t_out, y(:, 3), LineWidth=2)
legend('Noisy x', 'Actual x','Noisy y', 'Actual y','Noisy z', 'Actual z')
title('Magnetic Field Vector Components noisy vs actual (principal frame)')
xlabel('time [s]')
ylabel('Vector components')
%% check for eclipse condition
ind = find(eclipse_condition == 1);
disp((t_out(ind(end)) - t_out(ind(1)))/60)
disp((t_out(ind(end)) - t_out(ind(1)))/T_orbit)

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
%% Plotting Angular Momentum Vector
% Analytical Solution for axially symmetric satellite
%Angular Momentum
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
% plot(t_out, quat_out(:, 4), LineWidth=2)
% hold on;
% plot(t_out, quat_out(:, 1),  LineWidth=2)
% plot(t_out, quat_out(:, 2),  LineWidth=2)
% plot(t_out, quat_out(:, 3), LineWidth=2)
% xlabel('t [s]')
% ylabel('Quaternions')
% legend('q4', 'q1', 'q2', 'q3')
% title('Quaternions over 5 orbits')
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
% % 
% %%
% Plot Omega over time
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
% figure()
% hold on;
% plot(t_out, rad2deg(unwrap(euler_out(:, 1))), LineWidth=2)
% plot(t_out, rad2deg(unwrap(euler_out(:, 2))), LineWidth=2)
% plot(t_out, rad2deg(unwrap(euler_out(:, 3))), LineWidth=2)
% xlabel('t [s]')
% ylabel('Angle [deg]')
% legend('\phi', '\theta', '\psi')
% title('Euler angles over time (213 sequence)')
% 
% %Plot Omega over time
% figure()
% hold on;
% plot(t_out, omega_out_inertial(:, 1), LineWidth=2)
% plot(t_out, omega_out_inertial(:, 2), LineWidth=2)
% plot(t_out, omega_out_inertial(:, 3), LineWidth=2)
% xlabel('t [s]')
% ylabel('\omega [rad/s]')
% legend('X', 'Y', 'Z')
% title('Inertial angular velocity over time')


%% Problem 4
% Verify magnitude
% %Plot Torque over time
% figure()
% hold on;
% plot(t_out, torque_out(:, 1), LineWidth=2)
% plot(t_out, torque_out(:, 2),LineWidth=2)
% plot(t_out, torque_out(:, 3),LineWidth=2)
% plot(t_out, sqrt(torque_out(:, 1).^2 + torque_out(:, 2).^2 + torque_out(:, 3).^2), LineWidth=2)
% xlabel('t [s]')
% ylabel('M [Nm]')
% legend('M_x', 'M_y', 'M_z', 'M_{tot}')
% title('Torque over time (one orbit)')

%Plot Torque over time
% figure()
% hold on;
% plot(t_out, M_grav_out(:, 1), LineWidth=2)
% plot(t_out, M_grav_out(:, 2),LineWidth=2)
% plot(t_out, M_grav_out(:, 3),LineWidth=2)
% plot(t_out, sqrt(M_grav_out(:, 1).^2 + M_grav_out(:, 2).^2 + M_grav_out(:, 3).^2), LineWidth=2)
% yline(6.167e-5,'--k','LineWidth',2)
% xlabel('t [s]')
% ylabel('M [Nm]')
% legend('M_x', 'M_y', 'M_z', 'M_{tot}','M_{max}')
%title('Gravity gradient torque over time (one orbit)')

%% Plotting PSET 5

% % Verify magnetic torque
% figure()
% hold on;
% plot(t_out, M_mag_out(:, 1), LineWidth=2)
% plot(t_out, M_mag_out(:, 2),LineWidth=2)
% plot(t_out, M_mag_out(:, 3),LineWidth=2)
% plot(t_out, sqrt(M_mag_out(:, 1).^2 + M_mag_out(:, 2).^2 + M_mag_out(:, 3).^2), LineWidth=2)
% M_mag_max = 2*0.000154100322136014*(6378^3)*(3.08e-5)/((7.0049e3)^3)
% yline(M_mag_max,'--k','LineWidth',2)
% yline(-M_mag_max,'--k','LineWidth',2)
% xlabel('t [s]')
% ylabel('Magnetic Torque [Nm]')
% legend('M_x', 'M_y', 'M_z', 'M_{tot}','M_{max}')
% title('Magnetic torque over time (one orbit)')
% 
% % Verify SRP torque
% figure()
% hold on;
% plot(t_out, M_SRP_out(:, 1), LineWidth=2)
% plot(t_out, M_SRP_out(:, 2),LineWidth=2)
% plot(t_out, M_SRP_out(:, 3),LineWidth=2)
% plot(t_out, sqrt(M_SRP_out(:, 1).^2 + M_SRP_out(:, 2).^2 + M_SRP_out(:, 3).^2), LineWidth=2)
% M_SRP_max = (1358/(3e8))*((8660.25+9139.44+6*6908.59)/(100^2))*(1+0.85)*(15/100)
% % yline(M_SRP_max,'--k','LineWidth',2)
% xlabel('t [s]')
% ylabel('Solar Radiation Pressure Torque [Nm]')
% legend('M_x', 'M_y', 'M_z', 'M_{tot}','M_{max}')
% title('SRP torque over time (one orbit)')
% 
% % Verify Aero torque
% figure()
% hold on;
% plot(t_out, M_aero_out(:, 1), LineWidth=2)
% plot(t_out, M_aero_out(:, 2),LineWidth=2)
% plot(t_out, M_aero_out(:, 3),LineWidth=2)
% plot(t_out, sqrt(M_aero_out(:, 1).^2 + M_aero_out(:, 2).^2 + M_aero_out(:, 3).^2), LineWidth=2)
% M_aero_max = 0.5*7e-14*((7.5634e3)^2)*1.28*((8660.25+9139.44+6*6908.59)/(100^2))*(15/100)
% % yline(M_aero_max,'--k','LineWidth',2)
% xlabel('t [s]')
% ylabel('Aerodynamic Torque [Nm]')
% legend('M_x', 'M_y', 'M_z', 'M_{tot}','M_{max}')
% title('Aerodynamic torque over time (one orbit)')
% 
% % All torques
% figure()
% hold on;
% plot(t_out, sqrt(M_mag_out(:, 1).^2 + M_mag_out(:, 2).^2 + M_mag_out(:, 3).^2), LineWidth=2)
% plot(t_out, sqrt(M_SRP_out(:, 1).^2 + M_SRP_out(:, 2).^2 + M_SRP_out(:, 3).^2), LineWidth=2)
% plot(t_out, sqrt(M_aero_out(:, 1).^2 + M_aero_out(:, 2).^2 + M_aero_out(:, 3).^2), LineWidth=2)
% plot(t_out, sqrt(M_grav_out(:, 1).^2 + M_grav_out(:, 2).^2 + M_grav_out(:, 3).^2), LineWidth=2)
% xlabel('t [s]')
% ylabel('Torque [Nm]')
% legend('M_{mag}', 'M_{SRP}', 'M_{aero}', 'M_{grav}')
% title('All torques over time (one orbit)')
% 
% % Resultant torques
% figure()
% hold on;
% plot(t_out, M_ALL_out, LineWidth=2)
% xlabel('t [s]')
% ylabel('Torque [Nm]')
% legend('M_x', 'M_y', 'M_z')
% title('Resultant total torque over time (one orbit)')
% 
% 
% %% 
%% Plotting PSET6
%Plot Estimated VS actual attitude (EULER)
%% Plot
%Deterministic Attitude Determination Euler Angles over time
euler_est_det = DCMseries2eulerseries(DCM_estimated_det_out);
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_est_det(:, 1)), 'red', 'LineWidth',2)
plot(t_out, rad2deg(euler_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('\phi [deg]')
title('Deterministic Attitude Determination: Euler angles over time (213 sequence)')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_est_det(:, 2)), 'red',  'LineWidth',2)
plot(t_out, rad2deg(euler_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_est_det(:, 3)),'Color', 'red', 'LineWidth',2)
plot(t_out, rad2deg(euler_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('\psi [deg]')
legend('Deterministic Attitude Determination', 'True Euler Angles')
%% sign match quaternions
for i = 1:length(quat_out)
    if(sign(quat_out(i, 1)) ~= sign(quat_estimated_Q_out(i, 1)))
        quat_estimated_Q_out(i, :) = - quat_estimated_Q_out(i, :);
    end
end

%% Statistical
figure()
subplot(4,1,1)
hold on;
plot(t_out, (quat_estimated_Q_out(:, 1)), 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q1')
title('Quaternions estimated from statistical attitude determination')
subplot(4,1,2)
hold on;
plot(t_out, (quat_estimated_Q_out(:, 2)), 'red',  'LineWidth',2)
plot(t_out, (quat_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q2')
subplot(4,1,3)
hold on;
plot(t_out, (quat_estimated_Q_out(:, 3)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q3')
subplot(4,1,4)
hold on;
plot(t_out, (quat_estimated_Q_out(:, 4)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 4)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q4')
legend('Estimated Quaternions', 'True Quaternions')

%% Kinematic Estimation
figure()
subplot(4,1,1)
hold on;
plot(t_out, (quat_estimated_kin_out(:, 1)), 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q1')
title('Quaternions estimated from kinematic equations')
subplot(4,1,2)
hold on;
plot(t_out, (quat_estimated_kin_out(:, 2)), 'red',  'LineWidth',2)
plot(t_out, (quat_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q2')
subplot(4,1,3)
hold on;
plot(t_out, (quat_estimated_kin_out(:, 3)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q3')
subplot(4,1,4)
hold on;
plot(t_out, (quat_estimated_kin_out(:, 4)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 4)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q4')
legend('Estimated Quaternions', 'True Quaternions')


