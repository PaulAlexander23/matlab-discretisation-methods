classdef Domain
    properties
        x
        shape
        dimension
    end
    
    methods
        function obj = Domain(x)
            obj.x = x;
            obj.shape = cellfun(@length,obj.x)';
            obj.dimension = length(obj.shape);
        end
        
    end
    
end