%% return omega_dot (zero torque)
function [omegadot] = omega_dot(t, omega)
    %[R_princ,inertia_p] = inertia();
    inertia_p = [85.075, 0, 0; 0, 85.075, 0; 0, 0, 120.2515]; %axial
    %symmetric satellite
    M = zeros(3,1); %no torques
    %Euler eqns in principal axes
    omegadot = [ 1/inertia_p(1,1) * (M(1) - (inertia_p(3,3) - inertia_p(2,2)) * omega(2)*omega(3));
                 1/inertia_p(2,2) * (M(2) - (inertia_p(1,1) - inertia_p(3,3)) * omega(3)*omega(1));
                 1/inertia_p(3,3) * (M(3) - (inertia_p(2,2) - inertia_p(1,1)) * omega(1)*omega(2))];
end