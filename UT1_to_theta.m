function [theta] = UT1_to_theta(UT1)
    M = UT1(1);
    D = UT1(2);
    Y = UT1(3);
    H = UT1(4);
    if M <= 2
        y = Y-1;
        m = M+12;
    else
        y = Y;
        m = M;
    end
    if (y > 1582 || (y == 1582 && (m >= 10 || (m == 10 && D > 10))))
        B = (y/400) - (y/100) + (y/4);
    elseif (y == 1582 && m == 10 && D < 10 && D > 4)
        B = 0;
    else
        B = -2 + ((y + 4716)/4) - 1179;
    end
    MJD = 365*y - 679004 + floor(B) + floor(30.6001*(m+1)) + D;
    MJD = MJD + (H/24);

    d = MJD - 51544.5;
    GMST = 280.4606 + 360.9856473*d;
    GMST = GMST*pi/180;
    GMST = wrapTo2Pi(GMST);

    theta = rotz(-GMST*180/pi);
end