classdef Domain
    properties
        x
        shape
        dimension
    end
    
    methods
        function obj = Domain(x)
            obj.x = x;
            obj.dimension = length(obj.x);
            obj.shape = cellfun(@length,obj.x)';
            if length(obj.shape) == 1, obj.shape = [obj.shape, 1]; end
            
        end
        
    end
    
end