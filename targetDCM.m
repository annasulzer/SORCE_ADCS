function target_DCM = targetDCM(r_sun, R_princ)
    r_sun_norm = r_sun/norm(r_sun);

    %define that body x-axis will be in plane with body z-axis and r_sun
    %(perpendicular on r_sun)
    r_sun_projected = [r_sun_norm(1:2); 0];
    target_y_body = cross(r_sun_projected, r_sun_norm)./norm(cross(r_sun_projected, r_sun_norm));
    target_x_body = cross(target_y_body, r_sun_norm); %in plane

    %disp(dot(r_sun_norm, target_x_body)) %check that they are perependicular
    
    %construct target DCM in principal frame
    target_DCM = eye(3); %target DCM in principal frame
    target_DCM(1:3, 1) =  R_princ * target_x_body; %target attitude x-component principal frame
    target_DCM(1:3, 2) =  R_princ * target_y_body; %target attitude y-component principal frame
    target_DCM(1:3, 3) =  R_princ * r_sun_norm; %target attitude z-component principal frame
end