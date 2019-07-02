classdef FDDomain < Domain
    properties
        D
        diffDegrees
        accuracy
        BCs
    end
    
    methods
        function obj = FDDomain(x, diffDegrees, accuracy, BCs)
            obj = obj@Domain(x);
            if nargin >= 4
                obj.BCs = BCs;
            else
                obj.BCs = cell(obj.dimension, 1);
                obj.BCs(:) = {"Periodic"};
            end
            obj.diffDegrees = diffDegrees;
            obj.accuracy = accuracy;
            obj.D = initialiseDiffMatrices(obj, obj.diffDegrees);
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
        
        function D = initialiseDiffMatrices(obj, diffDegrees)
            D = cell(1, size(diffDegrees, 2));
            
            for j = 1:size(diffDegrees, 2)
                D{j} = 1;
                for dim = 1:obj.dimension
                    Dmat = obj.constructMatrix(obj.x{dim}, diffDegrees(dim, j), dim);
                    D{j} = kron(Dmat, D{j});
                end
            end
        end
        
        function index = findDiffCellIndex(obj, deg)
            [~, index] = ismember(deg(:, 1)', obj.diffDegrees', 'rows');
        end
        
    end
    
    methods(Access = private)
        function mat = constructMatrix(obj, x, degree, dim)
            [stencil, coefficients] = obj.centredScheme(degree, obj.accuracy);
            xLength = length(x);
            dx = x(2) - x(1);
            B = repmat(coefficients*dx.^(-degree), xLength, 1);
            mat = spdiags(B, stencil, xLength, xLength);
            mat = obj.addBoundaryConditions(mat, degree, dim, x, xLength);
        end
        
        function mat = addBoundaryConditions(obj, mat, degree, dim, x, xLength)
            if isstring(obj.BCs{dim})
                if obj.BCs{dim} == "Periodic"
                    [stencil, coefficients] = obj.centredScheme(degree, obj.accuracy);
                    dx = x(2) - x(1);
                    B = repmat(coefficients*dx.^(-degree), xLength, 1);
                    mat = mat + obj.spdiagsPeriodicBCs(B, stencil, xLength);
                end
            else
                mat = sparse(xLength,xLength);
            end
        end
        
        function [stencil, coefficients] = centredScheme(obj, degree, accuracy)
            stencilLength = ceil(degree/2) + ceil(accuracy/2) - 1;
            stencil = -stencilLength:stencilLength;
            ex = (0:2 * stencilLength)';
            S = repmat(stencil, stencilLength*2+1, 1).^ex;
            del = factorial(degree) * (ex == degree);
            coefficients = (S \ del)';
        end
        
        function mat = spdiagsPeriodicBCs(obj, B, stencil, size1)
            mat = spdiags(B, stencil+size1, size1, size1);
            mat = spdiags(B, stencil-size1, size1, size1) + mat;
        end
    end
end