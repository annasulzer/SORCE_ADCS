%% Quaternions
function quaternions_init = inital_quaternions(C)
    beta_4sq = (0.25*(1+trace(C)));
    beta_1sq = (0.25*(1+2*C(1,1)- trace(C)));
    beta_2sq = (0.25*(1+2*C(2,2)- trace(C)));
    beta_3sq = (0.25*(1+2*C(3,3)- trace(C)));
    
    [maximum, max_ind] = max([beta_1sq, beta_2sq, beta_3sq, beta_4sq]);
    
    switch max_ind
        case 1
            beta_1 = sqrt(beta_1sq);
            beta_4 = (C(2,3) - C(3,2))/(4*beta_1);
            beta_3 = (C(3,1) + C(1,3))/(4*beta_1);
            beta_2 = (C(1,2) + C(2,1))/(4*beta_1);
        case 2
            beta_2 = sqrt(beta_2sq);
            beta_4 = (C(3,1) - C(1,3))/(4*beta_2);
            beta_3 = (C(2,3) + C(3,2))/(4*beta_2);
            beta_1 = (C(1,2) + C(2,1))/(4*beta_2);
        case 3
            beta_3 = sqrt(beta_3sq);
            beta_4 = (C(1,2) - C(2,1))/(4*beta_3);
            beta_2 = (C(2,3) + C(3,2))/(4*beta_3);
            beta_1 = (C(3,1) + C(1,3))/(4*beta_3);
           
        case 4
            beta_4 = sqrt(beta_4sq);
            beta_3 = (C(1,2) - C(2,1))/(4*beta_4);
            beta_2 = (C(3,1) - C(1,3))/(4*beta_4);
            beta_1 = (C(2,3) - C(3,2))/(4*beta_4);   
    end
    
    quaternions_init = [beta_4; beta_1; beta_2; beta_3];
    quaternions_init =  quaternions_init/ norm(quaternions_init); %normalize
end

