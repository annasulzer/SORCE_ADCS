%% Returns RTN axis expressed in inertial frame. Input is the state vector in inertial
%  Anna Sulzer & Ethan Anzia
%  AA279C PSET3
function [R, T, N] = RTN_frame_inertial(state)
    n = length(state(:,1));
    R = zeros(n, 3);
    T = zeros(n, 3);
    N = zeros(n, 3);
    for i = 1:n
        R(i, :) = state(i, 1:3)/norm(state(i, 1:3));
        N(i, :) = cross(state(i, 1:3), state(i, 4:6))/norm(cross(state(i, 1:3),state(i, 4:6)));
        T(i, :) = cross(N(i, :), R(i, :))/norm(cross(N(i, :), R(i, :)));
    end
end