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


% Sun IC
UT1 = [1,25,2004,00];
JD_init = 2453029.5;
MJD_init = JD_init - 2400000.5;

D_init = JD_init - 2451545.0;
theta = UT1_to_theta(UT1);

% EKF
EKF_P0 = eye(7);
EKF_P0(1:3, 1:3) = 0.005.* EKF_P0(1:3, 1:3); %gyro variance
EKF_P0(4:7, 4:7) = 0.01.* EKF_P0(4:7, 4:7); %angle variance on sun

%Initial conditions
[state_ECI_init, T_orbit, n] = OrbitPropagation();
omega_init = [0; 0; 0.01];
DCM_initial = targetDCM([26321453.5527815,	-132781955.130633,	-57571626.5531097]', R_princ); %rinitial state sun


%Integration settings
eps = 1e-10;
absTol= 1e-10;
relTol = 1e-6;
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

% EKF
EKF_P_post_min = out.EKF_P_post_min.Data;
EKF_x_post_min = out.EKF_x_post_min.Data(:,:)';
EKF_P_post = out.EKF_P_post.Data;
EKF_x_post = out.EKF_x_post.Data(:,:)';

eclipse_condition = out.eclipse.Data();


% %% Plotting PSET6
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
% %% sign match quaternions
% for i = 1:length(quat_out)
%     if(sign(quat_out(i, 4)) ~= sign(quat_estimated_Q_out(i, 4)))
%         quat_estimated_Q_out(i, :) = - quat_estimated_Q_out(i, :);
%     end
% end
% 
% %% Statistical
% figure()
% subplot(4,1,1)
% hold on;
% plot(t_out, (quat_estimated_Q_out(:, 1)), 'red', 'LineWidth',2)
% plot(t_out, (quat_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q1')
% title('Quaternions estimated from statistical attitude determination')
% subplot(4,1,2)
% hold on;
% plot(t_out, (quat_estimated_Q_out(:, 2)), 'red',  'LineWidth',2)
% plot(t_out, (quat_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q2')
% subplot(4,1,3)
% hold on;
% plot(t_out, (quat_estimated_Q_out(:, 3)),'Color', 'red', 'LineWidth',2)
% plot(t_out, (quat_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q3')
% subplot(4,1,4)
% hold on;
% plot(t_out, (quat_estimated_Q_out(:, 4)),'Color', 'red', 'LineWidth',2)
% plot(t_out, (quat_out(:, 4)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q4')
% legend('Estimated Quaternions', 'True Quaternions')
% 
% %% Kinematic Estimation
% figure()
% subplot(4,1,1)
% hold on;
% plot(t_out, (quat_estimated_kin_out(:, 1)), 'red', 'LineWidth',2)
% plot(t_out, (quat_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q1')
% title('Quaternions estimated from kinematic equations')
% subplot(4,1,2)
% hold on;
% plot(t_out, (quat_estimated_kin_out(:, 2)), 'red',  'LineWidth',2)
% plot(t_out, (quat_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q2')
% subplot(4,1,3)
% hold on;
% plot(t_out, (quat_estimated_kin_out(:, 3)),'Color', 'red', 'LineWidth',2)
% plot(t_out, (quat_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q3')
% subplot(4,1,4)
% hold on;
% plot(t_out, (quat_estimated_kin_out(:, 4)),'Color', 'red', 'LineWidth',2)
% plot(t_out, (quat_out(:, 4)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
% xlabel('t [s]')
% ylabel('q4')
% legend('Estimated Quaternions', 'True Quaternions')


%% PSET7 Plotting
%% Deterministic
DCM_error_det_estimated = DCM_out;
for i = 1:length(DCM_out)
    DCM_error_det_estimated(:, :, i) = DCM_out(:,:,i) * DCM_estimated_det_out(:,:,i)';
end
euler_error_estimated_det = DCMseries2eulerseries(DCM_error_det_estimated);

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error_estimated_det(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Error in Euler Angles for Deterministic Attitude Determination')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error_estimated_det(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error_estimated_det(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')

%% Q method
DCM_error_Q_estimated = zeros(3,3,length(DCM_out));
DCM_estimated_Q = DCM_out;
for i = 1:length(DCM_out)
    DCM_estimated_Q(:,:,i) = quat2dcm([quat_estimated_Q_out(i, 4), quat_estimated_Q_out(i, 1:3)]);
    DCM_error_Q_estimated(:, :, i) = DCM_out(:,:,i) * DCM_estimated_Q(:,:,i)';
end
euler_error_estimated_Q = DCMseries2eulerseries(DCM_error_Q_estimated);

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error_estimated_Q(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Error in Euler Angles for Statistical Attitude Determination')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error_estimated_Q(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error_estimated_Q(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')

%% Kinematic
DCM_error_kin_estimated = zeros(3,3,length(DCM_out));
DCM_estimated_kin = DCM_out;
for i = 1:length(DCM_out)
    DCM_estimated_kin(:,:,i) = quat2dcm([quat_estimated_kin_out(i, 4), quat_estimated_kin_out(i, 1:3)]);
    DCM_error_kin_estimated(:, :, i) = DCM_out(:,:,i) * DCM_estimated_kin(:,:,i)';
end
euler_error_estimated_kin = DCMseries2eulerseries(DCM_error_kin_estimated);

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error_estimated_kin(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Error in Euler Angles for Kinematic Attitude Determination')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error_estimated_kin(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error_estimated_kin(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')

%% small Angles and mean/cov
V_det = var(rad2deg(euler_error_estimated_det));
mean_det = mean(rad2deg(euler_error_estimated_det));
V_Q = var(rad2deg(euler_error_estimated_Q));
mean_Q = mean(rad2deg(euler_error_estimated_Q));
V_kin = var(rad2deg(euler_error_estimated_kin));
mean_kin = mean(rad2deg(euler_error_estimated_kin));


DCM_error_Q_estimated_smallAngles = zeros(3,3,length(DCM_out));
for i = 1:length(DCM_out)
    DCM_error_Q_estimated_smallAngles(:, :, i) = DCM_error_Q_estimated(:,:,i) * small_angle_DCM(euler_error_estimated_Q(i, :))';
end
euler_error_estimated_Q_small_angle = DCMseries2eulerseries(DCM_error_Q_estimated_smallAngles);

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error_estimated_Q_small_angle(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Small Euler Angle Approximation Error for Q-Method')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error_estimated_Q_small_angle(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error_estimated_Q_small_angle(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')


DCM_error_det_estimated_smallAngles = zeros(3,3,length(DCM_out));
for i = 1:length(DCM_out)
    DCM_error_det_estimated_smallAngles(:, :, i) = DCM_error_det_estimated(:,:,i) * small_angle_DCM(euler_error_estimated_det(i, :))';
end
euler_error_estimated_det_small_angle = DCMseries2eulerseries(DCM_error_det_estimated_smallAngles);

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error_estimated_det_small_angle(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Small Euler Angle Approximation Error for Deterministic AD')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error_estimated_det_small_angle(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error_estimated_det_small_angle(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')


DCM_error_kin_estimated_smallAngles = zeros(3,3,length(DCM_out));
for i = 1:length(DCM_out)
    DCM_error_kin_estimated_smallAngles(:, :, i) = DCM_error_kin_estimated(:,:,i) * small_angle_DCM(euler_error_estimated_kin(i, :))';
end
euler_error_estimated_kin_small_angle = DCMseries2eulerseries(DCM_error_kin_estimated_smallAngles);

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error_estimated_kin_small_angle(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Small Euler Angle Approximation Error for kinematic AD')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error_estimated_kin_small_angle(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error_estimated_kin_small_angle(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')

%frobenius norms
frobenius_norm_smallAngles = zeros(3,length(DCM_out));
for i = 1:length(DCM_out)
    frobenius_norm_smallAngles(1, i) = norm(DCM_error_det_estimated(:,:,i) - small_angle_DCM(euler_error_estimated_det(i, :)), 'fro');
    frobenius_norm_smallAngles(2, i) = norm(DCM_error_Q_estimated(:,:,i) - small_angle_DCM(euler_error_estimated_Q(i, :)), 'fro');
    frobenius_norm_smallAngles(3, i) = norm(DCM_error_kin_estimated(:,:,i) - small_angle_DCM(euler_error_estimated_kin(i, :)), 'fro');
end

%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, frobenius_norm_smallAngles(1, :))
xlabel('t [s]')
ylabel('Frobenius Norm')
legend('Deterministic AD')
title('Frobenius Norm of the difference between small angle DCM and exact DCM')
subplot(3,1,2)
plot(t_out, frobenius_norm_smallAngles(2, :))
xlabel('t [s]')
ylabel('Frobenius Norm')
legend('Q-Method')
subplot(3,1,3)
plot(t_out, frobenius_norm_smallAngles(3, :))
xlabel('t [s]')
ylabel('Frobenius Norm')
legend('Kinematic AD')

%% EKF

figure()
subplot(4,1,1)
hold on;
plot(t_out, (EKF_x_post_min(:, 4)), 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q1')
title('Quaternions estimated from kinematic equations')
subplot(4,1,2)
hold on;
plot(t_out, (EKF_x_post_min(:, 5)), 'red',  'LineWidth',2)
plot(t_out, (quat_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q2')
subplot(4,1,3)
hold on;
plot(t_out, (EKF_x_post_min(:, 6)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q3')
subplot(4,1,4)
hold on;
plot(t_out, (EKF_x_post_min(:, 7)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 4)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q4')
legend('Estimated Quaternions', 'True Quaternions')

%errors in quat
figure()
subplot(4,1,1)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 1))-(EKF_x_post_min(20:end-2, 4)), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q1')
title('Quaternion Error between EKF estimated and true (with perturbations')
subplot(4,1,2)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 2))-(EKF_x_post_min(20:end-2, 5)), 'red',  'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q2')
subplot(4,1,3)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 3))-(EKF_x_post_min(20:end-2, 6)),'Color', 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q3')
subplot(4,1,4)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 4))-(EKF_x_post_min(20:end-2, 7)),'Color', 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q4')



%% errors in omega
figure()
subplot(3,1,1)
hold on;
plot(t_out(20:end-2), (omega_out(20:end-2, 1))-(EKF_x_post_min(20:end-2, 1)), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta \omega_x')
title('Angular Velocity Error EKF VS. True (with perturbations')
subplot(3,1,2)
hold on;
plot(t_out(20:end-2), (omega_out(20:end-2, 2))-(EKF_x_post_min(20:end-2, 2)), 'red',  'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta \omega_y')
subplot(3,1,3)
hold on;
plot(t_out(20:end-2), (omega_out(20:end-2, 3))-(EKF_x_post_min(20:end-2, 3)),'Color', 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta \omega_z')


%%
% Plot over time
figure()
subplot(3,1,1)
plot(t_out, EKF_x_post_min(:, 1),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, omega_out(:, 1), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('\omega_x [rad/s]')
title('Angular Velocity over time from the EKF VS true state')
subplot(3,1,2)
plot(t_out, EKF_x_post_min(:, 2),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, omega_out(:, 2), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)

xlabel('t [s]')
ylabel('\omega_y [rad/s]')
subplot(3,1,3)
plot(t_out, EKF_x_post_min(:, 3),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, omega_out(:, 3), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)

xlabel('t [s]')
ylabel('\omega_z [rad/s]')
legend('Estimated Omega', 'True Omega')
%% Uncertainty plots
std_wx = squeeze(sqrt(EKF_P_post_min(1,1,:)));
std_wy = squeeze(sqrt(EKF_P_post_min(2,2,:)));
std_wz = squeeze(sqrt(EKF_P_post_min(3,3,:)));
std_q1 = squeeze(sqrt(EKF_P_post_min(4,4,:)));
std_q2 = squeeze(sqrt(EKF_P_post_min(5,5,:)));
std_q3 = squeeze(sqrt(EKF_P_post_min(6,6,:)));
std_q4 = squeeze(sqrt(EKF_P_post_min(7,7,:)));


figure()
subplot(3,1,1)
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,1) + std_wx; flipud(EKF_x_post_min(:,1) - std_wx)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);
xlabel('t [s]')
ylabel('\omega_x [rad/s]')
title('Angular Velocity over time with the standard deviations from the EKF')
subplot(3,1,2)
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,2) + std_wy; flipud(EKF_x_post_min(:,2) - std_wy)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);

xlabel('t [s]')
ylabel('\omega_y [rad/s]')
subplot(3,1,3)
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,3) + std_wz; flipud(EKF_x_post_min(:,3) - std_wz)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);

xlabel('t [s]')
ylabel('\omega_z [rad/s]')
legend('Standard Deviation', 'Mean Omega')

%errors in quat
figure()
subplot(4,1,1)
hold on;
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,4) + std_q1; flipud(EKF_x_post_min(:,4) - std_q1)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);
xlabel('t [s]')
ylabel('q1')
title('Quaternion over time with the standard deviations from the EKF')
subplot(4,1,2)
hold on;
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,5) + std_q2; flipud(EKF_x_post_min(:,5) - std_q2)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);
xlabel('t [s]')
ylabel('q2')
subplot(4,1,3)
hold on;
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,6) + std_q3; flipud(EKF_x_post_min(:,6) - std_q3)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);
xlabel('t [s]')
ylabel('q3')
subplot(4,1,4)
hold on;
fill([t_out; flipud(t_out)],[EKF_x_post_min(:,7) + std_q4; flipud(EKF_x_post_min(:,7) - std_q4)],'g', 'FaceAlpha', 0.15, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, EKF_x_post_min(:,1), 'b', 'LineWidth', 2);
xlabel('t [s]')
ylabel('q4')
legend('Standard Deviation', 'Mean Quaternion')

%%
% Plot over time
figure()
subplot(4,1,1)
plot(t_out, EKF_x_post_min(:, 1),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, omega_out(:, 1), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('\q1 [rad/s]')
title('Angular Velocity over time from the EKF VS true state')
subplot(4,1,2)
plot(t_out, EKF_x_post_min(:, 2),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, omega_out(:, 2), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)

xlabel('t [s]')
ylabel('\q2 [rad/s]')
subplot(4,1,3)
plot(t_out, EKF_x_post_min(:, 3),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, omega_out(:, 3), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)

xlabel('t [s]')
ylabel('\q3 [rad/s]')
subplot(4,1,3)
plot(t_out, EKF_x_post_min(:, 3),'Color', 'red', 'LineWidth',2)
hold on;
plot(t_out, quat_out(:, 3), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)

xlabel('t [s]')
ylabel('\q4 [rad/s]')
legend('Estimated Quat', 'True Quat')
%% calculate small angle DCM for 213 sequence
function DCM = small_angle_DCM(angles)
    ay = angles(1);%phi
    ax = angles(2); %theta
    az = angles(3);
    DCM = [1, az, -ay; -az, 1, ax; ay, -ax, 1]; 
end

