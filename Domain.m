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

        function Y = reshapeToVector(obj, y)
            ycell = obj.extractSurfacesAsCell(y);

            Ycell = cellfun(@(y) reshape(y, [prod(obj.shape), 1]), ycell, ...
                'UniformOutput', false);

            Y = cell2mat(Ycell);
        end

        function y = reshapeToDomain(obj, Y)
            Ycell = obj.extractVectorAsCell(Y);

            ycell = cellfun(@(Y) reshape(Y, obj.shape), Ycell, ...
                'UniformOutput', false);

            y = cell2mat(ycell);
        end

        function fq = interp(obj, x, f, method)
            if nargin < 4; method = 'linear'; end

            fq = interpn(x{:}, f, obj.x{:}, method);
        end

        function ycell = extractSurfacesAsCell(obj, y)
            ycell = obj.extractShapeAsCell(y, obj.shape);
        end

        function Ycell = extractVectorAsCell(obj, Y)
            Ycell = obj.extractShapeAsCell(Y, [prod(obj.shape), 1]);
        end
    end

    methods (Access = private)
        function shape = calculateShape(~, x)
            shape = reshape(cellfun(@length, x), 1, []);
            if length(shape) == 1
                shape = [shape, 1];
            end
        end

        function y = reshapeIndependentVariables(obj,x)
            xShapes = ones(obj.dimension) + eye(obj.dimension) .* (obj.shape - 1);
            y = cell(obj.dimension, 1);
            for dim = 1:obj.dimension
                y{dim} = reshape(x{dim}, xShapes(dim,:));
            end
        end

        function ycell = extractShapeAsCell(~, y, shape)
            shapeMultiplier = size(y)./shape;
            ycell = mat2cell(y, ...
                repmat(shape(1), 1, shapeMultiplier(1)), ...
                repmat(shape(2), 1, shapeMultiplier(2)));
        end
    end
end
