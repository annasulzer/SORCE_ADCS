%% Euler Equations
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2
%% 
clear; clc; close all;
%%
%Initial conditions
 omega_init = [0.1; 0.15; 0.08];
% omega_init = [0.25; 0.001; 0.001];
% omega_init = [0.001; 0.25; 0.001];
% omega_init = [0.001; 0.001; 0.25];

%[R_princ,inertia_p] = inertia(); %inertia in principal axes
inertia_p = [85.075, 0, 0; 0, 85.075, 0; 0, 0, 120.2515];

time = 2*pi/norm(omega_init);

%Integrate
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
tstart = 0; tint = 0.01; tend = 100*time;
[t_out, omega_out] = ode113(@omega_dot,[tstart:tint:tend]', omega_init, options);
L_out = (inertia_p * omega_out')';

% Plot over time
figure()
subplot(3,1,1)
plot(t_out, omega_out(:, 1))
xlabel('t [s]')
ylabel('\omega_x [rad/s]')
title('Angular Velocity over time')
subplot(3,1,2)
plot(t_out, omega_out(:, 2))
xlabel('t [s]')
ylabel('\omega_y [rad/s]')
subplot(3,1,3)
plot(t_out, omega_out(:, 3))
xlabel('t [s]')
ylabel('\omega_z [rad/s]')


%saveas(gcf, 'ang_vel.pdf');
%% Ellipsoids

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


% Analytical Solutions Ellipsoid
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

%Angular Momentum
L_analyt=zeros(length(t_out), 3);
L_analyt_init = inertia_p * omega_init;
for i = 1:length(t_out)
    L_analyt(i, 3) = L_analyt_init(3); %omega_z 
    L_analyt(i, 1) = norm(L_analyt_init(1:2)) * cos(lambda*t_out(i) + theta_zero); %omega_x = om_xy * cos(lambda*t)
    L_analyt(i, 2) = norm(L_analyt_init(1:2)) * sin(lambda*t_out(i) + theta_zero); %omega_x = om_xy * sin(lambda*t)
end

%% Plot
%Ellipsoid
figure()
surface(x_enE, y_enE, z_enE, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);
hold on;
surface(x_momE, y_momE, z_momE, 'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha',.2, 'EdgeAlpha',.4);

%Polhode
plot3(omega_out(:,1), omega_out(:,2), omega_out(:, 3), 'Color', 'black', 'LineWidth', 2)
plot3(omega_analyt(:,1), omega_analyt(:,2), omega_analyt(:, 3), 'Color', 'green', 'LineWidth', 2)

axis equal;
grid on;
xlabel('\omega_x [1/s^2]')
ylabel('\omega_y [1/s^2]')
zlabel('\omega_z [1/s^2]')
legend('Energy Ellipsoid', 'Momentum Ellipsoid', 'Polhode Numerical', 'Polhode Analytical')
view(3)
title('Analytical vs Numerical Polhode for axial-symmetric satellite')

%Errors btw analytical and numerical
errors = omega_out - omega_analyt;

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


%Plot Ang Vel over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, omega_out(:, 1), 'black')
plot(t_out, omega_analyt(:, 1), '--','Color','green')
xlabel('t [s]')
ylabel('\omega_x [rad/s]')
title('Numerical solution and analytical angular velocity over time')
subplot(3,1,2)
hold on;
plot(t_out, omega_out(:, 2), 'black')
plot(t_out, omega_analyt(:, 2), '--','Color','green')
xlabel('t [s]')
ylabel('\omega_y [rad/s]')
subplot(3,1,3)
hold on;
plot(t_out, omega_out(:, 3), 'black')
plot(t_out, omega_analyt(:, 3),'--', 'Color','green')
xlabel('t [s]')
ylabel('\omega_z [rad/s]')



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

%% Quaternions
C = [0.892539, 0.157379, -0.422618;
    -0.275451, 0.932257, -0.234570;
    0.357073, 0.325773, 0.875426];
beta_4sq = (0.25*(1+trace(C)));
beta_1sq = (0.25*(1+2*C(1,1)- trace(C)));
beta_2sq = (0.25*(1+2*C(2,2)- trace(C)));
beta_3sq = (0.25*(1+2*C(3,3)- trace(C)));

[maximum, max_ind] = max([beta_1sq, beta_2sq, beta_3sq, beta_4sq]);

switch max_ind
    case 1
        beta_1 = sqrt(beta_1sq);
        beta_4 = (C(2,3) - C(3,2))/(4*beta_1);
        beta_3 = (C(3,1) + C(1,3))/(4*beta_1);
        beta_2 = (C(1,2) + C(2,1))/(4*beta_1);
    case 2
        beta_2 = sqrt(beta_2sq);
        beta_4 = (C(3,1) - C(1,3))/(4*beta_2);
        beta_3 = (C(2,3) + C(3,2))/(4*beta_2);
        beta_1 = (C(1,2) + C(2,1))/(4*beta_2);
    case 3
        beta_3 = sqrt(beta_3sq);
        beta_4 = (C(1,2) - C(2,1))/(4*beta_3);
        beta_2 = (C(2,3) + C(3,2))/(4*beta_3);
        beta_1 = (C(3,1) + C(1,3))/(4*beta_3);
       
    case 4
        beta_4 = sqrt(beta_4sq);
        beta_3 = (C(1,2) - C(2,1))/(4*beta_4);
        beta_2 = (C(3,1) - C(1,3))/(4*beta_4);
        beta_1 = (C(2,3) - C(3,2))/(4*beta_4);   
end

quaternions_init = [beta_4; beta_1; beta_2; beta_3];
quaternions_init =  quaternions_init/ norm(quaternions_init); %normalize

%% Kinematic Equation of motion
function [quat_dot] = quaternions_dot(quaternion, omega)
    Sigma = [0, omega(3), -omega(2), omega(1);
            -omega(3), 0, omega(1), omega(2);
            omega(2), -omega(1), 0, omega(3);
            -omega(1), -omega(2), -omega(3), 0];
    quat_dot = 0.5 * Sigma * quaternion;
end

%% return omega_dot (zero torque)
function [omegadot] = omega_dot(t, omega)
    %[R_princ,inertia_p] = inertia();
    inertia_p = [85.075, 0, 0; 0, 85.075, 0; 0, 0, 120.2515];
    M = zeros(3,1); %no torques
    %Euler eqns in principal axes
    omegadot = [ 1/inertia_p(1,1) * (M(1) - (inertia_p(3,3) - inertia_p(2,2)) * omega(2)*omega(3));
                 1/inertia_p(2,2) * (M(2) - (inertia_p(1,1) - inertia_p(3,3)) * omega(3)*omega(1));
                 1/inertia_p(3,3) * (M(3) - (inertia_p(2,2) - inertia_p(1,1)) * omega(1)*omega(2))];
end