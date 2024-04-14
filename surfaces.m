% SORCE SADC
% External Surfaces Discritization 
% Ethan Anzia and Anna Sulzer

function [surface] = getsurface(surfacenum)
    if surfacenum == 1
        surface.bary = [-1.29,-0.19,100.22];
        surface.size = 8660.25;
        surface.unit = [0,0,1];
    end
    if surfacenum == 2
        surface.bary = [-1.29,-0.19,-60.08];
        surface.size = 8660.25;
        surface.unit = [0,0,-1];
    end
    if surfacenum == 3
        surface.bary = [48.71,-0.19,20.07];
        surface.size = 9139.44;
        surface.unit = [1,0,0];
    end
    if surfacenum == 4
        surface.bary = [108.54,-0.19,-58.08];
        surface.size = 6908.59;
        surface.unit = [0,0,1];
    end
    if surfacenum == 5
        surface.bary = [108.54,-0.19,-60.08];
        surface.size = 6908.59;
        surface.unit = [0,0,-1];
    end
    if surfacenum == 6
        surface.bary = [23.71,43.11,20.07];
        surface.size = 9139.44;
        surface.unit = [0.5,0.866,0];
    end
    if surfacenum == 7
        surface.bary = [53.62,94.92,-58.08];
        surface.size = 6908.59;
        surface.unit = [0,0,1];
    end
    if surfacenum == 8
        surface.bary = [53.62,94.92,-60.08];
        surface.size = 6908.59;
        surface.unit = [0,0,-1];
    end
    if surfacenum == 9
        surface.bary = [-26.29,43.11,20.07];
        surface.size = 9139.44;
        surface.unit = [-0.5,0.866,0];
    end
    if surfacenum == 10
        surface.bary = [-56.20,94.92,-58.08];
        surface.size = 6908.59;
        surface.unit = [0,0,1];
    end
    if surfacenum == 11
        surface.bary = [-56.20,94.92,-60.08];
        surface.size = 6908.59;
        surface.unit = [0,0,-1];
    end
    if surfacenum == 12
        surface.bary = [-51.29,-0.19,20.07];
        surface.size = 9139.44;
        surface.unit = [-1,0,0];
    end
    if surfacenum == 13
        surface.bary = [-111.11,-0.19,-58.08];
        surface.size = 6908.59;
        surface.unit = [0,0,1];
    end
    if surfacenum == 14
        surface.bary = [-111.11,-0.19,-60.08];
        surface.size = 6908.59;
        surface.unit = [0,0,-1];
    end
    if surfacenum == 15
        surface.bary = [-26.29,-43.49,20.07];
        surface.size = 9139.44;
        surface.unit = [-0.5,-0.866,0];
    end
    if surfacenum == 16
        surface.bary = [-56.20,-95.30,-58.08];
        surface.size = 6908.59;
        surface.unit = [0,0,1];
    end
    if surfacenum == 17
        surface.bary = [-56.20,-95.30,-60.08];
        surface.size = 6908.59;
        surface.unit = [0,0,-1];
    end
    if surfacenum == 18
        surface.bary = [23.71,-43.49,20.07];
        surface.size = 9139.44;
        surface.unit = [0.5,-0.866,0];
    end
    if surfacenum == 19
        surface.bary = [53.62,-95.30,-58.08];
        surface.size = 6908.59;
        surface.unit = [0,0,1];
    end
    if surfacenum == 20
        surface.bary = [53.62,-95.30,-60.08];
        surface.size = 6908.59;
        surface.unit = [0,0,-1];
    end
end



