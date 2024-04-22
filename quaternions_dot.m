%% Kinematic Equation of motion
function [quat_dot] = quaternions_dot(quaternion, omega)
    Sigma = [0, omega(3), -omega(2), omega(1);
            -omega(3), 0, omega(1), omega(2);
            omega(2), -omega(1), 0, omega(3);
            -omega(1), -omega(2), -omega(3), 0];
    quat_dot = 0.5 * Sigma * quaternion;
end
