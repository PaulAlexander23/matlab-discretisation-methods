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
            obj.shape = reshape(cellfun(@length,obj.x), 1, []);
            if length(obj.shape) == 1, obj.shape = [obj.shape, 1]; end
        end
        
        function Y = reshapeToVector(obj, y)
            Y = reshape(y, prod(obj.shape), []);
        end
        
        function y = reshapeToDomain(obj, Y)
            if obj.dimension == 1
                y = reshape(Y, [obj.shape(1), numel(Y)/obj.shape(1)]);
            else
                y = reshape(Y, [obj.shape, numel(Y)/prod(obj.shape)]);
            end
        end
    end
    
end