%% Orbit Propagator
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2
%%
clear; clc;
% Initial conditions
mu = 398600.435436; %Earth centered orbit
rE = 6378;
e = 0.0026520;
om = 337.7127;%deg
Om = 11.4426;%deg
i = 39.9975; %deg
n = 14.80808427* 2*pi/86400; % mean motion in rad per sec
a = (mu/(n^2))^(1/3);

%Perifocal
r = a*(1-e); %at periapsis
v = sqrt(2*mu/r - mu/a); %circular orbit speed
state_PQW_init = [r, 0, 0, 0, v, 0, 0]; %at periapsis

%ECI
R_PQW_IJK = rotz(Om)*rotx(i)*rotz(om); %rotation matrix from PQW to IJK, swapped sign bc active vs passive rotation
state_ECI_init = [R_PQW_IJK*state_PQW_init(1:3)',R_PQW_IJK*state_PQW_init(4:6)'] ;
T_orbit = 2*pi*sqrt(a^3/mu);

%Orbit Propagator ode113
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
tstart = 0; tint = 0.01; tend = 10*T_orbit;
[t_out, state_out] = ode113(@state_dot,[tstart:tint:tend]', state_ECI_init, options);

%% Plot
figure()
plot3(state_out(:, 1), state_out(:, 2), state_out(:, 3), 'red', 'LineWidth',2)

%Earth
opts.FaceAlpha = 0.6;
opts.Units = 'km';
opts.RotAngle = -127.1871; 
opts.RefPlane = 'ecliptic';
planet3D('Earth',opts);
axis equal;

view(3);
grid on;
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Orbit around a Earth in ECI')
hold off;

%% return equations of motion
function [statedot] = state_dot(t, state)
    mu = 398600.435436; %Earth centered orbit
    statedot = zeros(6, 1);
    statedot(1:3) = state(4:6);
    statedot(4:6) = -mu * state(1:3)/(norm(state(1:3))^3);
end