%% CAD with Principal Axes for PSET 2
clear; clc;

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






