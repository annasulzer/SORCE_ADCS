function [R, T, N] = RTN_frame_inertial(state)
    n = length(state(:,1));
    R = zeros(3, n);
    T = zeros(3, n);
    N = zeros(3, n);
    for i = 1:T
        R(i) = state(i, 1:3)/norm(state(i, 1:3));
        N(i) = cross(state(i, 1:3),state(i, 1:3))/norm(cross(state(i, 1:3),state(i, 1:3)));
        T(i) = cross(N, R)/norm(cross(N,R));
    end
    R = R';
    T = T';
    N = N';
end