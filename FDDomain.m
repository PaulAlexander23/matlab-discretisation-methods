classdef FDDomain < Domain
    properties
        D
        diffDegrees
        accuracy
        direction
    end

    methods
        function obj = FDDomain(x, diffDegrees, accuracy, direction)
            if nargin < 4, direction = "central"; end
            
            obj = obj@Domain(x);
            
            obj.direction = direction;
            obj.diffDegrees = diffDegrees;
            obj.accuracy = accuracy;
            obj.D = initialiseDiffMatrices(obj, obj.diffDegrees, obj.accuracy);
            
        end

        function dy = diff(obj, y, degree)
            
            dy = reshapeToVector(obj, y);

            Dx = obj.diffMat(degree);

            largeD = kron(speye(size(dy,1)/size(Dx,1)), Dx);
            
            dy = largeD * dy;
            
            dy = reshapeToDomain(obj, dy);
        end
        
        function D = diffMat(obj, degree)
            index = findDiffCellIndex(obj,degree);
            D = obj.D{index};
        end
        
        function w = multiply(obj, u, v, powers)
            if nargin < 4
                powers = [1, 1];
            end
            w = u.^powers(1) .* v.^powers(2);
        end
    end

    methods(Access = private)
        function D = initialiseDiffMatrices(obj, diffDegrees, accuracy)
            D = cell(1, size(diffDegrees, 2));

            for j = 1:size(diffDegrees, 2)
                D{j} = 1;
                for k = 1:obj.dimension
                    Dmat = constructMatrix(obj, obj.x{k}, diffDegrees(k, j), accuracy);
                    D{j} = kron(Dmat, D{j});
                end
            end

            function mat = constructMatrix(obj, x, degree, accuracy)
                if obj.direction == "central"
                    [stencil, coefficients] = centredScheme(degree, accuracy);
                elseif obj.direction == "forward"
                    [stencil, coefficients] = forwardScheme(degree, accuracy);
                else
                    error('Unimplemented direction: %s, use central or forward',...
                        obj.direction);
                end
                xLength = length(x);
                dx = x(2) - x(1);
                B = repmat(coefficients*dx.^(-degree), xLength, 1);
                mat = spdiagsPeriodicBCs(B, stencil, xLength);
            end

            function [stencil, coefficients] = centredScheme(degree, accuracy)
                stencilLength = ceil(degree/2) + ceil(accuracy/2) - 1;
                stencil = -stencilLength:stencilLength;
                ex = (0:2 * stencilLength)';
                S = repmat(stencil, stencilLength*2+1, 1).^ex;
                del = factorial(degree) * (ex == degree);
                coefficients = (S \ del)';
            end
            
            function [stencil, coefficients] = forwardScheme(degree, accuracy)
                stencilLength = degree + accuracy;
                stencil = 0:stencilLength - 1;
                ex = (0:stencilLength-1)';
                S = repmat(stencil, stencilLength, 1).^ex;
                del = factorial(degree) * (ex == degree);
                coefficients = (S \ del)';
            end

            function mat = spdiagsPeriodicBCs(B, stencil, size1)
                mat = spdiags(B, stencil, size1, size1);
                mat = spdiags(B, stencil+size1, size1, size1) + mat;
                mat = spdiags(B, stencil-size1, size1, size1) + mat;
            end
        end

        function index = findDiffCellIndex(obj, deg)
            [~, index] = ismember(deg(:, 1)', obj.diffDegrees', 'rows');
        end
    end
end
