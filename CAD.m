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
% Run EulerEquations 

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

t = 0:1:300;
title(sprintf('Angular Momentum Vector Over Time\nTime: %0.2f sec', t(1)), 'Interpreter','latex','FontSize',16);
filename = 'Ang_Momentum_Vector.gif';

quiver3(161.35+0.19,0+60.08,169.65+1.29,L_out(1,1),L_out(1,2),L_out(1,3),LineWidth=3,ShowArrowHead='on',Color='red',MaxHeadSize=1)

for k = 1:5:length(t)
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









