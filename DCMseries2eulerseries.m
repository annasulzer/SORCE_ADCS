function euler_series = DCMseries2eulerseries(DCM)
    euler_series = zeros(length(DCM(1,1,:)), 3);
    for i = 1:length(DCM(1, 1, :))
        %get euler angles 213
        euler_series(i, 1) = atan2(-DCM(1, 3, i), DCM(3, 3, i));
        euler_series(i, 2) = asin(DCM(2, 3, i));
        euler_series(i, 3) = atan2(-DCM(2, 1, i), DCM(2, 2, i));
    end
end