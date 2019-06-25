classdef FDDomain < Domain
    properties
        D
        diffDegrees
        accuracy
    end

    methods
        function obj = FDDomain(x, diffDegrees, accuracy)
            obj = obj@Domain(x);

            obj.diffDegrees = diffDegrees;
            obj.accuracy = accuracy;
            obj.D = initialiseDiffMatrices(obj, obj.diffDegrees, obj.accuracy);
        end

        function dy = diff(obj, y, degree)
            dy = reshapeToVector(obj, y);

            index = findDiffCellIndex(obj, degree);

            dy = obj.D{index} * dy;

            dy = reshapeToDomain(obj, dy);
        end

        function D = getDiffMatrix(obj, deg)
            D = obj.D{findDiffCellIndex(obj, deg)};
        end
    end

    methods(Access = private)
        function D = initialiseDiffMatrices(obj, diffDegrees, accuracy)
            D = cell(1, size(diffDegrees, 2));

            for j = 1:size(diffDegrees, 2)
                D{j} = 1;
                for k = 1:obj.dimension
                    Dmat = constructMatrix(obj.x{k}, diffDegrees(k, j), accuracy);
                    D{j} = kron(Dmat, D{j});
                end
            end

            function mat = constructMatrix(x, degree, accuracy)
                [stencil, coefficients] = centredScheme(degree, accuracy);
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