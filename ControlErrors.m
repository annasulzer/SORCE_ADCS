%% Control Errors
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET6

%Control Reference: Sun Vecotr aligned with Body-z-axis
alpha = 45; %azimuth from x-axis sun
beta = 10; %elevation sun
r_sun_norm = [cosd(alpha)*cosd(beta); sind(alpha)*cosd(beta); sind(beta)]; %target for body z-axis

%define that body x-axis will be in plane with body z-axis and r_sun
%(perpendicular on r_sun)
target_x_body = [-cosd(alpha)*sind(beta); -sind(alpha)*sind(beta); cosd(beta)];
disp(dot(r_sun_norm, target_x_body)) %check that they are perependicular

%construct target DCM in principal frame
target_DCM = eye(3); %target DCM in principal frame
target_DCM(1:3, 1) =  R_princ * target_x_body; %target attitude x-component principal frame
target_DCM(1:3, 2) =  R_princ * cross(r_sun_norm, target_x_body); %target attitude y-component principal frame
target_DCM(1:3, 3) =  R_princ * r_sun_norm; %target attitude z-component principal frame

%%%%%%%%%%%%%%%%%Simulate with Euler Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Errors
DCM_error = DCM_out;
euler_error = zeros(length(DCM_out(1,1,:)), 3);
for i = 1:length(DCM_out)
    DCM_error(:, :, i) = target_DCM * DCM_out(:, :, i)';
    %get euler angles 213
    euler_error(i, 1) = atan2(-DCM_error(1, 3, i), DCM_error(3, 3, i));
    euler_error(i, 2) = asin(DCM_error(2, 3, i));
    euler_error(i, 3) = atan2(-DCM_error(2, 1, i), DCM_error(2, 2, i));
end

%Validate
alignment_error = zeros(length(DCM_out(1,1,:)), 3);
for i = 1:length(DCM_out)
    DCM_alignment_error = DCM_error(:, :, i) * DCM_out(:, :, i);
    %get euler angles 213
    alignment_error(i, 1) = atan2(-DCM_alignment_error(1, 3), DCM_alignment_error(3, 3));
    alignment_error(i, 2) = asin(DCM_alignment_error(2, 3));
    alignment_error(i, 3) = atan2(-DCM_alignment_error(2, 1), DCM_alignment_error(2, 2));
end

%% Plot
%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(euler_error(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Euler angles over time (213 sequence)')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(euler_error(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(euler_error(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')


%% Plot alignment error
%Euler Angles Errors over time
figure()
subplot(3,1,1)
hold on;
plot(t_out, rad2deg(alignment_error(:, 1)), 'blue')
xlabel('t [s]')
ylabel('\Delta\phi [deg]')
title('Alignment Error over time (Euler 213 sequence)')
subplot(3,1,2)
hold on;
plot(t_out, rad2deg(alignment_error(:, 2)), 'blue')
xlabel('t [s]')
ylabel('\Delta\theta [deg]')
subplot(3,1,3)
hold on;
plot(t_out, rad2deg(alignment_error(:, 3)), 'blue')
xlabel('t [s]')
ylabel('\Delta\psi [deg]')
