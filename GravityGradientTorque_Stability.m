%% Gravity Gradient Torque Stability Calculation
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET5
clear; clc; close all;
[R_princ,inertia_p] = inertia(); %inertia in principal axes
%inertia_p = diag([120, 85, 60]); %unstable in all
%inertia_p = diag([100, 85, 110]); %unstable pitch
%inertia_p = diag([60, 150, 120]); %unstable roll yaw


%stability coefficients
KN = (inertia_p(2,2) - inertia_p(1,1))/ inertia_p(3,3);
KR = (inertia_p(3,3) - inertia_p(2,2))/ inertia_p(1,1);%Kt
KT = (inertia_p(3,3) - inertia_p(1,1))/ inertia_p(2,2);%Kr

%specify grid
kR = -1:0.001:1 ; 
kT = -1:0.001:1 ; 
[kT_m, kR_m] = meshgrid(kR, kT);

%stability conditions
cond1 = kT_m>kR_m;
cond2 = kR_m.*kT_m >0;
cond3 = abs(1+3*kT_m+kR_m.*kT_m) >= 4*sqrt(kR_m.*kT_m) ;

%calculate where condition fulfilled
region1_false = nan(size(kR_m));
region2_false = nan(size(kR_m));
region3_false = nan(size(kR_m));


region1_false(~cond1) = 1;
region2_false(~cond2) = 1;
region3_false(~cond3) = 1;


% Plot
figure()
hold on
imagesc(kR, kT,3*double(region3_false),'AlphaData', 1.0); %blue
imagesc(kR, kT,3*double(region2_false),'AlphaData', 0.5); %blue
imagesc(kR, kT,2*double(region1_false),'AlphaData', 0.45);%yellow, pitch unstable
imagesc(kR, kT,double((region2_false.*region1_false)),'AlphaData', .25); %green
imagesc(kR, kT,double((region3_false.*region1_false)),'AlphaData', .25); %green

colormap([1 1 1; 18/255 135/255 49/255; 250/255 189/255 35/255; 18/255 94/255 135/255])%white, green, yellow, blue
set(gca, 'CLim', [0 3]);

%satellite
plot(KT, KR, 'x', 'MarkerSize',5, 'Color','black', 'LineWidth',2)

%for legend
fill([nan nan nan], [nan nan nan], [18/255 94/255 135/255], 'FaceAlpha', 0.3)
fill([nan nan nan], [nan nan nan], [250/255 189/255 35/255], 'FaceAlpha', 0.3)
fill([nan nan nan], [nan nan nan], [18/255 135/255 49/255], 'FaceAlpha', 0.3)

%show axis
plot([0 0], ylim, 'k',  'HandleVisibility','off'); % Vertical line at x=0
plot(xlim, [0 0], 'k',  'HandleVisibility','off'); % Horizontal line at y=0


xlabel('K_T')
ylabel('K_R')
axis equal
grid on
xlim([-1,1])
ylim([-1,1])
legend('Satellite', 'Unstable Roll and Yaw', 'Unstable Pitch', 'Unstable Roll, Pitch and Yaw', '', '')

