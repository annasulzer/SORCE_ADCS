%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2 and PSET3
%  For Polhode and Ellipsoid Plotting run EllipsoidAndPolhode.m
%% 
clear; clc; close all;
%%
%Initial conditions
 omega_init = [0.1; 0.15; 0.08];

 % Quaternions
DCM = [0.892539, 0.157379, -0.422618;
     -0.275451, 0.932257, -0.234570;
     0.357073, 0.325773, 0.875426];

%[R_princ,inertia_p] = inertia(); %inertia in principal axes
inertia_p = [85.075, 0, 0; 0, 85.075, 0; 0, 0, 120.2515]; %axial symmetric

%Integration settings
absTol= 1e-10;
relTol = 1e-6;
time = 2*pi/norm(omega_init);
tstart = 0; tend = 100*time;


%%
out = sim("main");
omega_out = out.omega.Data(:,:)';
t_out = out.tout;
L_out = (inertia_p * omega_out')';
quat_out = out.quaternions.Data(:, :)';
%% Plotting
% Plot over time
figure()
subplot(3,1,1)
plot(t_out, omega_out(:, 1))
xlabel('t [s]')
ylabel('\omega_x [rad/s]')
title('Angular Velocity over time from the numerical integration')
subplot(3,1,2)
plot(t_out, omega_out(:, 2))
xlabel('t [s]')
ylabel('\omega_y [rad/s]')
subplot(3,1,3)
plot(t_out, omega_out(:, 3))
xlabel('t [s]')
ylabel('\omega_z [rad/s]')

%saveas(gcf, 'ang_vel.pdf');

%% Analytical Solution for axially symmetric satellite
%Angular Velocity
omega_analyt=zeros(length(t_out), 3);
lambda = omega_init(3)*(inertia_p(3,3) - inertia_p(1,1))/inertia_p(1,1);
theta_zero = atan2(omega_init(2), omega_init(1));
for i = 1:length(t_out)
    omega_analyt(i, 3) = omega_init(3); %omega_z 
    omega_analyt(i, 1) = norm(omega_init(1:2)) * cos(lambda*t_out(i) + theta_zero); %omega_x = om_xy * cos(lambda*t)
    omega_analyt(i, 2) = norm(omega_init(1:2)) * sin(lambda*t_out(i) + theta_zero); %omega_x = om_xy * sin(lambda*t)
end
%Errors btw analytical and numerical
errors = omega_out - omega_analyt;

%Angular Momentum
L_analyt=zeros(length(t_out), 3);
L_analyt_init = inertia_p * omega_init;
for i = 1:length(t_out)
    L_analyt(i, 3) = L_analyt_init(3); %omega_z 
    L_analyt(i, 1) = norm(L_analyt_init(1:2)) * cos(lambda*t_out(i) + theta_zero); %omega_x = om_xy * cos(lambda*t)
    L_analyt(i, 2) = norm(L_analyt_init(1:2)) * sin(lambda*t_out(i) + theta_zero); %omega_x = om_xy * sin(lambda*t)
end

% Plot over time
figure()
subplot(3,1,1)
plot(t_out, omega_analyt(:, 1))
xlabel('t [s]')
ylabel('\omega_x [rad/s]')
title('Angular Velocity over time from the analytical solution')
subplot(3,1,2)
plot(t_out, omega_analyt(:, 2))
xlabel('t [s]')
ylabel('\omega_y [rad/s]')
subplot(3,1,3)
plot(t_out, omega_analyt(:, 3))
xlabel('t [s]')
ylabel('\omega_z [rad/s]')

%Plot Error over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, errors(:, 1), 'red')
xlabel('t [s]')
ylabel('\Delta\omega_x [rad/s]')
title('Errors between numerical solution and analytical solution over time')
subplot(3,1,2)
hold on;
plot(t_out, errors(:, 2), 'red')
xlabel('t [s]')
ylabel('\Delta\omega_y [rad/s]')
subplot(3,1,3)
hold on;
plot(t_out, errors(:, 3), 'red')
xlabel('t [s]')
ylabel('\Delta\omega_z [rad/s]')


%Plot Angular Momentum
figure()
subplot(3,1,1)
hold on;
plot(t_out, L_out(:, 1), 'black')
plot(t_out, L_analyt(:, 1),'--', 'Color','green')
xlabel('t [s]')
ylabel('L_x [kg*m^2/s]')
title('Angular momentum vector of numerical solution and analytical solution over time')
subplot(3,1,2)
hold on;
plot(t_out, L_out(:, 2), 'black')
plot(t_out, L_analyt(:, 2), '--','Color','green')
xlabel('t [s]')
ylabel('L_y [kg*m^2/s]')
subplot(3,1,3)
hold on;

plot(t_out, L_out(:, 3), 'black')
plot(t_out, L_analyt(:, 3),'--', 'Color','green')
xlabel('t [s]')
ylabel('L_z [kg*m^2/s]')



