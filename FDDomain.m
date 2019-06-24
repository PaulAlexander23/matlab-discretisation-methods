classdef FDDomain
    properties
        x
        D
        diffDegrees
    end
    
    methods
        function obj = FDDomain(x, degree, accuracy)
            obj.x = x;
            obj.D = init_fd(x, degree, accuracy);
            obj.diffDegrees = degree;
        end
        
        function dy = diff(obj, y, deg)
            dy = diff_fd(obj.x, y, deg, obj.D, obj.diffDegrees);
        end
    end
    
end