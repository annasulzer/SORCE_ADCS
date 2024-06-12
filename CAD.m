%% CAD with Principal Axes for PSET 2

filename = 'SORCE_forSTL.STL';

gm = importGeometry(filename);

figure
pdegplot(gm,"FaceAlpha",0.1)

delete(findobj(gca,'type','Text')); 
delete(findobj(gca,'type','Quiver'));

[principal_directions,principal_I_tensor] = inertia();

vector_rot_y = principal_directions*[50,0,0]';
vector_rot_z = principal_directions*[0,50,0]';
vector_rot_x = principal_directions*[0,0,50]';

hold on
quiver3(161.35,0,169.65,50,0,0,LineWidth=2,ShowArrowHead='on',Color='red')
quiver3(161.35,0,169.65,0,50,0,LineWidth=2,ShowArrowHead='on',Color='green')
quiver3(161.35,0,169.65,0,0,50,LineWidth=2,ShowArrowHead='on',Color='blue')

quiver3(161.35+0.19,0+60.08,169.65+1.29,50,0,0,LineWidth=2,ShowArrowHead='on',Color='red')
quiver3(161.35+0.19,0+60.08,169.65+1.29,0,50,0,LineWidth=2,ShowArrowHead='on',Color='green')
quiver3(161.35+0.19,0+60.08,169.65+1.29,0,0,50,LineWidth=2,ShowArrowHead='on',Color='blue')

quiver3(161.35+0.19,0+60.08,169.65+1.29,vector_rot_y(1),vector_rot_y(2),vector_rot_y(3),LineWidth=2,ShowArrowHead='on',Color='yellow')
quiver3(161.35+0.19,0+60.08,169.65+1.29,vector_rot_z(1),vector_rot_z(2),vector_rot_z(3),LineWidth=2,ShowArrowHead='on',Color='cyan')
quiver3(161.35+0.19,0+60.08,169.65+1.29,vector_rot_x(1),vector_rot_x(2),vector_rot_x(3),LineWidth=2,ShowArrowHead='on',Color='magenta')

legend('','','Y_{CAD}','Z_{CAD}','X_{CAD}','Y_{BODY}','Z_{BODY}','X_{BODY}','Z_{PRINCIPAL}','Y_{PRINCIPAL}(max)','X_{PRINCIPAL}(min)')

%% Angular Momentum Over Time
% Run EulerEquations First

filename = 'SORCE_forSTL.STL';
gm = importGeometry(filename);
figure
set(gcf,'position',[10,10,2000,800])
pdegplot(gm,"FaceAlpha",0.05)
hold on
delete(findobj(gca,'type','Text')); 
delete(findobj(gca,'type','Quiver'));

[principal_directions,principal_I_tensor] = inertia();

vector_rot_y = principal_directions*[50,0,0]';
vector_rot_z = principal_directions*[0,50,0]';
vector_rot_x = principal_directions*[0,0,50]';

t = t_out;
title(sprintf('Angular Momentum Vector Over Time\nTime: %0.2f sec', t(1)), 'Interpreter','latex','FontSize',16);
filename = 'Ang_Momentum_Vector.gif';

quiver3(161.35+0.19,0+60.08,169.65+1.29,L_out(1,1),L_out(1,2),L_out(1,3),LineWidth=3,ShowArrowHead='on',Color='red',MaxHeadSize=1)

for k = 1:1:length(t)
    view(135,15)
    delete(findobj(gca,'type','Text')); 
    delete(findobj(gca,'type','Quiver'));
    hold on
    quiver3(161.35+0.19,0+60.08,169.65+1.29,vector_rot_y(1),vector_rot_y(2),vector_rot_y(3),LineWidth=2,ShowArrowHead='on',Color='yellow')
    quiver3(161.35+0.19,0+60.08,169.65+1.29,vector_rot_z(1),vector_rot_z(2),vector_rot_z(3),LineWidth=2,ShowArrowHead='on',Color='cyan')
    quiver3(161.35+0.19,0+60.08,169.65+1.29,vector_rot_x(1),vector_rot_x(2),vector_rot_x(3),LineWidth=2,ShowArrowHead='on',Color='magenta')
    quiver3(161.35+0.19,0+60.08,169.65+1.29,2*L_out(k,1),2*L_out(k,2),2*L_out(k,3),LineWidth=3,ShowArrowHead='on',Color='red')
    quiver3(161.35+0.19,0+60.08,169.65+1.29,200*omega_out(k,1),200*omega_out(k,2),200*omega_out(k,3),LineWidth=3,ShowArrowHead='on',Color='blue')
    title(sprintf('Angular Momentum Vector Over Time\nTime: %0.2f sec', t(k)),...
    'Interpreter','Latex');
    legend('','','Z_{PRINCIPAL}','Y_{PRINCIPAL}(max)','X_{PRINCIPAL}(min)','Angular Momentum Vector (Scaled x2)','Angular Velocity Vector (Scaled x200)')
    pause(0.0000001)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end
end

%% Full Orbit Visualization
% Run main_script First! 
% Create Figure
figure; hold on 
title(sprintf('Orbit Propagation with RTN, Body, and Principal Axes\nTime: %0.2f sec', t_out(1)), 'Interpreter','latex','FontSize',16);
xlabel('x [km]','Interpreter','latex','FontSize',16)
ylabel('y [km]','Interpreter','latex','FontSize',16)
zlabel('z [km]','Interpreter','latex','FontSize',16)
view(3);
grid on;
set(gcf,'position',[10,10,2000,800]);

% Plot Earth
opts.FaceAlpha = 0.3;
opts.Units = 'km';
opts.RotAngle = -127.1871; 
opts.RefPlane = 'ecliptic';
planet3D('Earth',opts);
axis equal;

% Reassign from State Vector 
x_position = state_out(:,1);
y_position = state_out(:,2);
z_position = state_out(:,3);

% Create File 
filename = 'Orbit_Prop.gif';

% Plotting the first iteration
plot3(x_position,y_position,z_position,'Color','none');
xlim([-7000,7000])
ylim([-7000,7000])
zlim([-7000,7000])
p = plot3(x_position(1),y_position(1),z_position(1),'k','LineWidth',3);
m = scatter3(x_position(1),y_position(1),z_position(1),'filled','k');

% Iterating through the length of the time array
for k = 1:1:length(t_out)
    % Updating the line
    p.XData = x_position(1:k);
    p.YData = y_position(1:k);
    p.ZData = z_position(1:k);
    % Updating the point
    m.XData = x_position(k); 
    m.YData = y_position(k);
    m.ZData = z_position(k);
    % Delete Previous Quivers
    delete(findobj(gca,'type','Text')); 
    delete(findobj(gca,'type','Quiver'));
    hold on

    % RTN 
    quiver3(x_position(k),y_position(k),z_position(k),2000*R(k,1),2000*R(k,2),2000*R(k,3),LineWidth=3,ShowArrowHead='on',Color='#FFAC33')
    quiver3(x_position(k),y_position(k),z_position(k),2000*T(k,1),2000*T(k,2),2000*T(k,3),LineWidth=3,ShowArrowHead='on',Color='red')
    quiver3(x_position(k),y_position(k),z_position(k),2000*N(k,1),2000*N(k,2),2000*N(k,3),LineWidth=3,ShowArrowHead='on',Color='#7B3F00')

    % Principal
    quiver3(x_position(k),y_position(k),z_position(k),2000*Xp(k,1),2000*Xp(k,2),2000*Xp(k,3),LineWidth=3,ShowArrowHead='on',Color='magenta')
    quiver3(x_position(k),y_position(k),z_position(k),2000*Yp(k,1),2000*Yp(k,2),2000*Yp(k,3),LineWidth=3,ShowArrowHead='on',Color='cyan')
    quiver3(x_position(k),y_position(k),z_position(k),2000*Zp(k,1),2000*Zp(k,2),2000*Zp(k,3),LineWidth=3,ShowArrowHead='on',Color='yellow')

    % Body
    quiver3(x_position(k),y_position(k),z_position(k),2000*Xb(k,1),2000*Xb(k,2),2000*Xb(k,3),LineWidth=3,ShowArrowHead='on',Color='blue')
    quiver3(x_position(k),y_position(k),z_position(k),2000*Yb(k,1),2000*Yb(k,2),2000*Yb(k,3),LineWidth=3,ShowArrowHead='on',Color='red')
    quiver3(x_position(k),y_position(k),z_position(k),2000*Zb(k,1),2000*Zb(k,2),2000*Zb(k,3),LineWidth=3,ShowArrowHead='on',Color='green')

    % Legend
    legend('','','Orbit Path','','R','T','N','X_{principal}','Y_{principal}','Z_{principal}','X_{body}','Y_{body}','Z_{body}')

    % Updating the title
    title(sprintf('SORCE Orbit with RTN, Body, and Principal Axes\nTime: %0.2f sec', t_out(k)),...
    'Interpreter','Latex');
    % Delay
    pause(0.000000001)
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end
end

% function [statedot] = state_dot(t, state)
%     mu = 398600.435436; %Earth centered orbit
%     statedot = zeros(6, 1);
%     statedot(1:3) = state(4:6);
%     statedot(4:6) = -mu * state(1:3)/(norm(state(1:3))^3);
% end




