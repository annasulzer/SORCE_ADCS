function [Xp, Yp, Zp, Xb, Yb, Zb] = principal_frame_inert(R, R_princ)
    Xp = R(:,1,:)
    Yp = R(:,2,:)
    Zp = R(:,3,:)
    Xb = R_princ' * Xp;
    Yb = R_princ' * Yp;
    Zb = R_princ' * Zp;
end