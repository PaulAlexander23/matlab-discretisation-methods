classdef PSDomain < Domain
    properties
        
    end
    
    methods
        function obj = PSDomain(x)
            obj = obj@Domain(x);
        end
        
        function dy = diff(obj, y, deg)
            dy = diff_ps(obj.x, y, deg);
        end
    end
    
end