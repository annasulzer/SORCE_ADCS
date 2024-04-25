%% Returns Principal and body axis expressed in inertial frame. Input is 
%% Rotationmatrix from principal to inertial over time and R from body to principal
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET3

function [Xp, Yp, Zp, Xb, Yb, Zb] = principal_frame_inert(Rot, R_princ)
    n = length(Rot(1,1,:));
    Xp = zeros(n, 3);
    Yp = zeros(n, 3);
    Zp = zeros(n, 3);
    Xb = zeros(n, 3);
    Yb = zeros(n, 3);
    Zb = zeros(n, 3);
    for i = 1:n
        R = Rot(:, :, i);
        Xp(i, :) = R'*[1,0,0]';
        Yp(i, :) = R'*[0,1,0]';
        Zp(i, :) = R'*[0,0,1]';
        Xb(i, :) = R'* R_princ' * [1, 0, 0]';
        Yb(i, :) = R'* R_princ' * [0, 1, 0]';
        Zb(i, :) = R'* R_princ' * [0, 0, 1]';
    end
end