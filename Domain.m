classdef Domain
    properties
        x
        shape
        dimension
    end
    
    methods
        function obj = Domain(x)
            obj.dimension = length(x);
            obj.shape = obj.calculateShape(x);
            obj.x = obj.reshapeIndependentVariables(x);
        end
        
        function shape = calculateShape(~, x)
            shape = reshape(cellfun(@length, x), 1, []);
            if length(shape) == 1, shape = [shape, 1]; end
        end

        function y = reshapeIndependentVariables(obj,x)
            y = cell(obj.dimension, 1);
            xShapes = ones(obj.dimension) + eye(obj.dimension) .* (obj.shape - 1);
            for dim = 1:obj.dimension
                y{dim} = reshape(x{dim},xShapes(dim,:));
            end
        end
        
        function [Y, shapeMultiplier] = reshapeToVector(obj, y)
            
            shapeMultiplier = size(y)./obj.shape;
            
            if shapeMultiplier(1) == 1
                Y = reshape(y, prod(obj.shape), []);
            else
                ycell = mat2cell(y, ...
                    repmat(obj.shape(1),1,shapeMultiplier(1)), ...
                    repmat(obj.shape(2),1,shapeMultiplier(2)));
            
                Ycell = cellfun(@(y)reshape(y, prod(obj.shape),[]), ycell, ...
                    'UniformOutput', false);
            
                Y = cell2mat(Ycell);
            end
        end
        
        function [y, shapeMultiplier] = reshapeToDomain(obj, Y)
            
            shapeMultiplier = size(Y)./[prod(obj.shape),1];
            horizontalMultiplier = ones(1,length(obj.shape));
            horizontalMultiplier(2) = shapeMultiplier(2);
            
            y = squeeze(num2cell(reshape(Y, [obj.shape .* horizontalMultiplier,...
                shapeMultiplier(1)]),[1,2]));
            y = cat(1, y{:});
        end
        
        function fq = interp(obj, x, f)
            fq = interpn(x{:}, f, obj.x{:}, 'spine');
        end
    end
end
