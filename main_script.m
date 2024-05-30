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
EKF_P0(1:3, 1:3) = 0.005.* 10*EKF_P0(1:3, 1:3); %gyro variance *10
EKF_P0(4:7, 4:7) = 0.01.* 10* EKF_P0(4:7, 4:7); %angle variance on sun * 10

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
out = sim("mainEthan");
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
EKF_x_post = out.EKF_x_post.Data(:,:);
prefit_res = out.prefit_res.Data;
postfit_res = out.postfit_res.Data;

eclipse_condition = out.eclipse.Data();


%% PSET8 Plotting
%% Residuals Prefit
%Stats
mean_prefit_sun = mean(squeeze(prefit_res(:, 1, :)), 2)
cov_prefit_sun = cov(squeeze(prefit_res(:, 1, :))')
mean_prefit_mag = mean(squeeze(prefit_res(:, 2, :)), 2)
cov_prefit_mag = cov(squeeze(prefit_res(:, 2, :))')
mean_prefit_star = mean(squeeze(prefit_res(:, 3, :)), 2)
cov_prefit_star = cov(squeeze(prefit_res(:, 3, :))')

figure()
subplot(3,1,1)
hold on;
plot(t_out, squeeze(prefit_res(1, 1, :)), 'LineWidth',1)
ylim([-0.05,0.05])
xlabel('t [s]')
ylabel('z_x')
title('Prefit residuals Sun Sensor')
subplot(3,1,2)
hold on;
plot(t_out, squeeze(prefit_res(2, 1, :)), 'LineWidth',1)
ylim([-0.05,0.05])
xlabel('t [s]')
ylabel('z_y')
subplot(3,1,3)
hold on;
plot(t_out, squeeze(prefit_res(3, 1, :)), 'LineWidth',1)
ylim([-0.05,0.05])
xlabel('t [s]')
ylabel('z_z')

%Magnet
figure()
subplot(3,1,1)
hold on;
plot(t_out, squeeze(prefit_res(1, 2, :)), 'LineWidth',1)
ylim([-0.0000001, 0.0000001])
xlabel('t [s]')
ylabel('z_x [T]')
title('Prefit residuals Magnetometer')
subplot(3,1,2)
hold on;
plot(t_out, squeeze(prefit_res(2, 2, :)), 'LineWidth',1)
ylim([-0.0000001, 0.0000001])
xlabel('t [s]')
ylabel('z_y [T]')
subplot(3,1,3)
hold on;
plot(t_out, squeeze(prefit_res(3, 2, :)), 'LineWidth',1)
ylim([-0.0000001, 0.0000001])
xlabel('t [s]')
ylabel('z_z [T]')

%Star
figure()
subplot(3,1,1)
hold on;
plot(t_out, squeeze(prefit_res(1, 3, :)), 'LineWidth',1)
ylim([-0.001, 0.001])
xlabel('t [s]')
ylabel('z_x')
title('Prefit residuals Star Tracker (mean)')
subplot(3,1,2)
hold on;
plot(t_out, squeeze(prefit_res(2, 3, :)), 'LineWidth',1)
ylim([-0.001, 0.001])
xlabel('t [s]')
ylabel('z_y')
subplot(3,1,3)
hold on;
plot(t_out, squeeze(prefit_res(3, 3, :)), 'LineWidth',1)
ylim([-0.001, 0.001])
xlabel('t [s]')
ylabel('z_z')

%% Residuals Post fit
%Stats:
mean_postfit_sun = mean(squeeze(postfit_res(:, 1, :)), 2)
cov_postfit_sun = cov(squeeze(postfit_res(:, 1, :))')
mean_postfit_mag = mean(squeeze(postfit_res(:, 2, :)), 2)
cov_postfit_mag = cov(squeeze(postfit_res(:, 2, :))')
mean_postfit_star = mean(squeeze(postfit_res(:, 3, :)), 2)
cov_postfit_star = cov(squeeze(postfit_res(:, 3, :))')

% Residuals
figure()
subplot(3,1,1)
hold on;
plot(t_out, squeeze(postfit_res(1, 1, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_x')
title('Postfit residuals Sun Sensor')
subplot(3,1,2)
hold on;
plot(t_out, squeeze(postfit_res(2, 1, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_y')
subplot(3,1,3)
hold on;
plot(t_out, squeeze(postfit_res(3, 1, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_z')

%Magnet
figure()
subplot(3,1,1)
hold on;
plot(t_out, squeeze(postfit_res(1, 2, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_x [T]')
title('Postfit residuals Magnetometer')
subplot(3,1,2)
hold on;
plot(t_out, squeeze(postfit_res(2, 2, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_y [T]')
subplot(3,1,3)
hold on;
plot(t_out, squeeze(postfit_res(3, 2, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_z [T]')

%Star
figure()
subplot(3,1,1)
hold on;
plot(t_out, squeeze(postfit_res(1, 3, :)), 'LineWidth',1)
ylim([-5*10^(-4), 5*10^(-4)])
xlabel('t [s]')
ylabel('z_x')
title('Postfit residuals Star Tracker (mean)')
subplot(3,1,2)
hold on;
plot(t_out, squeeze(postfit_res(2, 3, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_y')
subplot(3,1,3)
hold on;
plot(t_out, squeeze(postfit_res(3, 3, :)), 'LineWidth',1)
xlabel('t [s]')
ylabel('z_z')
%% EKF
figure()
subplot(4,1,1)
hold on;
plot(t_out, (EKF_x_post(:, 4)), 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 1)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q1')
title('Quaternions estimated from kinematic equations')
subplot(4,1,2)
hold on;
plot(t_out, (EKF_x_post(:, 5)), 'red',  'LineWidth',2)
plot(t_out, (quat_out(:, 2)), 'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q2')
subplot(4,1,3)
hold on;
plot(t_out, (EKF_x_post(:, 6)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 3)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q3')
subplot(4,1,4)
hold on;
plot(t_out, (EKF_x_post(:, 7)),'Color', 'red', 'LineWidth',2)
plot(t_out, (quat_out(:, 4)),'Linestyle', '--', 'Color','blue', 'LineWidth',2)
xlabel('t [s]')
ylabel('q4')
legend('Estimated Quaternions', 'True Quaternions')

%% Calculate errors in quat
DCM_estimated = zeros(3,3,length(DCM_out));
errors_quat = zeros(size(EKF_x_post(:, 4:7)));
for i = 4:length(DCM_out)
    q4 = EKF_x_post(i, 7);
    q1 = EKF_x_post(i, 4);
    q2 = EKF_x_post(i, 5);
    q3 = EKF_x_post(i, 6);
     
    DCM_estimated = quat2dcm([EKF_x_post(i, 7), EKF_x_post(i, 4:6)]);
    DCM_error_quat = DCM_out(:,:,i) * DCM_estimated';
    errors = dcm2quat(DCM_error_quat);
    errors = errors./norm(errors);
    
    errors_quat(i, :) = [errors(2:4), errors(1)];
end


%% Uncertainty plots
std_wx = squeeze(sqrt(EKF_P_post(1,1,:)));
std_wy = squeeze(sqrt(EKF_P_post(2,2,:)));
std_wz = squeeze(sqrt(EKF_P_post(3,3,:)));
std_q1 = 3*squeeze(sqrt(EKF_P_post(4,4,:)));
std_q2 = 3*squeeze(sqrt(EKF_P_post(5,5,:)));
std_q3 = 3*squeeze(sqrt(EKF_P_post(6,6,:)));
std_q4 = 3*squeeze(sqrt(EKF_P_post(7,7,:)));

%Statistics
mean_errors_quat = mean(errors_quat)
cov_errors_quat = cov(errors_quat)
errors_om = (omega_out)-(EKF_x_post(:, 1:3));
mean_errors_om = mean(errors_om)
cov_errors_om = cov(errors_om)


figure()
subplot(3,1,1)
fill([t_out; flipud(t_out)],[errors_om(:,1) + std_wx; flipud(errors_om(:,1) - std_wx)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');

hold on;
plot(t_out, errors_om(:,1), 'b', 'LineWidth', 2);
ylim([-0.1, 0.1])

xlabel('t [s]')
ylabel('\Delta \omega_x [rad/s]')
title('Errors in Angular Velocity over time with the standard deviations from the EKF')
subplot(3,1,2)
fill([t_out; flipud(t_out)],[errors_om(:,2) + std_wy; flipud(errors_om(:,2) - std_wy)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, errors_om(:,1), 'b', 'LineWidth', 2);
ylim([-0.1, 0.1])

xlabel('t [s]')
ylabel('\Delta \omega_y [rad/s]')
subplot(3,1,3)
fill([t_out; flipud(t_out)],[errors_om(:,3) + std_wz; flipud(errors_om(:,3) - std_wz)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, errors_om(:,1), 'b', 'LineWidth', 2);
ylim([-0.1, 0.1])

xlabel('t [s]')
ylabel('\Delta \omega_z [rad/s]')
legend('Standard Deviation', 'True Error Omega')

%% errors in quat
figure()
subplot(4,1,1)
hold on;
fill([t_out; flipud(t_out)],[mean_errors_quat(:,1) + std_q1; flipud(mean_errors_quat(:,1) - std_q1)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, errors_quat(:,1), 'b', 'LineWidth', 1);
ylim([-0.0005, 0.0005])
xlabel('t [s]')
ylabel('\Delta q1')
title('Errors in Quaternions over time with the 3\sigma deviations from the EKF')

subplot(4,1,2)
hold on;
fill([t_out; flipud(t_out)],[mean_errors_quat(:,2) + std_q2; flipud(mean_errors_quat(:,2) - std_q2)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, errors_quat(:,2), 'b', 'LineWidth', 1);
ylim([-0.0005, 0.0005])
xlabel('t [s]')
ylabel('\Delta q2')

subplot(4,1,3)
hold on;
fill([t_out; flipud(t_out)],[mean_errors_quat(:,3) + std_q3; flipud(mean_errors_quat(:,3) - std_q3)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, errors_quat(:,3), 'b', 'LineWidth', 1);
ylim([0.0001, 0.0009])
xlabel('t [s]')
ylabel('\Delta q3')

subplot(4,1,4)
hold on;
fill([t_out; flipud(t_out)],[mean_errors_quat(:,4) + std_q4; flipud(mean_errors_quat(:,4) - std_q4)],'g', 'FaceAlpha', 1, 'FaceColor', '#8fce00', 'LineStyle', 'none');
hold on;
plot(t_out, errors_quat(:,4), 'b', 'LineWidth', 1);
ylim([0.9995, 1.0005])
xlabel('t [s]')
ylabel('\Delta q4')
legend('3\sigma Deviation', 'True Error Quaternion')

%% errors in quat
figure()
subplot(4,1,1)
hold on;
plot(t_out, (quat_out)-(EKF_x_post(:, 4)), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q1')
title('Quaternion Error between EKF estimated and true (with perturbations')
subplot(4,1,2)
hold on;
plot(t_out, (quat_out)-(EKF_x_post(:, 4)), 'red', 'LineWidth',1)
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


figure()
subplot(4,1,1)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 1))-(EKF_x_post(20:end-2, 4)), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q1')
title('Quaternion Error between EKF estimated and true (with perturbations')
subplot(4,1,2)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 2))-(EKF_x_post(20:end-2, 5)), 'red',  'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q2')
subplot(4,1,3)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 3))-(EKF_x_post(20:end-2, 6)),'Color', 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q3')
subplot(4,1,4)
hold on;
plot(t_out(20:end-2), (quat_out(20:end-2, 4))-(EKF_x_post(20:end-2, 7)),'Color', 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q4')


%%
figure()
subplot(4,1,1)
hold on;
plot(t_out, EKF_x_post(:, 4) - EKF_x_post_min(:, 4), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q1')
title('Quaternion Error between EKF prior and posterior')
subplot(4,1,2)
hold on;
plot(t_out, EKF_x_post(:, 5) - EKF_x_post_min(:, 5), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q2')
subplot(4,1,3)
hold on;
plot(t_out, EKF_x_post(:, 6) - EKF_x_post_min(:, 6), 'red', 'LineWidth',1)
xlabel('t [s]')
ylabel('\Delta q3')
subplot(4,1,4)
hold on;
plot(t_out, EKF_x_post(:, 7) - EKF_x_post_min(:, 7), 'red', 'LineWidth',1)
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
std_wx = squeeze(sqrt(EKF_P_post(1,1,:)));
std_wy = squeeze(sqrt(EKF_P_post(2,2,:)));
std_wz = squeeze(sqrt(EKF_P_post(3,3,:)));
std_q1 = squeeze(sqrt(EKF_P_post(4,4,:)));
std_q2 = squeeze(sqrt(EKF_P_post(5,5,:)));
std_q3 = squeeze(sqrt(EKF_P_post(6,6,:)));
std_q4 = squeeze(sqrt(EKF_P_post(7,7,:)));


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

