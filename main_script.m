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

%Initial conditions
[state_ECI_init, T_orbit, n] = OrbitPropagation();
omega_init = [0; 0; 0.0];

% Sun IC
UT1 = [1,25,2004,00];
JD_init = 2453029.5;
MJD_init = JD_init - 2400000.5;

D_init = JD_init - 2451545.0;
theta = UT1_to_theta(UT1);


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


eclipse_condition = out.eclipse.Data();
%% check that noise 
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
% Rot = DCM_out;
% [Xp, Yp, Zp, Xb, Yb, Zb] = principal_body_frame_inertial(Rot, R_princ);
% [R, T, N] = RTN_frame_inertial(state_out);


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


%% small Angles and mean/cov
V_det = var(rad2deg(euler_error_estimated_det));
mean_det = mean(rad2deg(euler_error_estimated_det));
V_Q = var(rad2deg(euler_error_estimated_Q));
mean_Q = mean(rad2deg(euler_error_estimated_Q));


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


function DCM = small_angle_DCM(angles)
    ay = angles(1);%phi
    ax = angles(2); %theta
    az = angles(3);
    DCM = [1, az, -ay; -az, 1, ax; ay, -ax, 1]; 
end

