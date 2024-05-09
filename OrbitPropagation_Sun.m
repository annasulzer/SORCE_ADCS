%% Orbit Propagator initializer for Sun
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET2


%orbit elements of sun around Earth
%(https://www.karhukoti.com/SunKep#google_vignette)
function [t0_MJD, a, Om, e, om, i, M0, n] = initialConditions_sun()
    t0_MJD = 24001.00000000000; %MJD 2024
    a = 149597870;
    Om = 0;
    e =  0.0167133;
    om = 282.7685;
    i = 23.4406;
    M0 = deg2rad(357.6205);
    n = 0.002737778522 * 2*pi/86400; %or sidereal day??
end

