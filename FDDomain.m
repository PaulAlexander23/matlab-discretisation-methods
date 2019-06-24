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
            obj.D = getDiffMatrices(obj, obj.diffDegrees, obj.accuracy);
        end
        
        function dy = diff(obj, y, deg)
            dy = reshapeToVector(obj, y);
            
            index = findDiffCellIndex(obj, deg);
            
            dy = obj.D{index} * dy;
            
            dy = reshapeToShape(obj, dy);
            
            dy = {dy};
        end
    end
    
    methods (Access=private)
        function D = getDiffMatrices(obj, diffDegrees, accuracy)
            
            D = cell(1,size(diffDegrees,2));
            
            for j = 1:size(diffDegrees,2)
                D{j} = 1;
                for k = 1:obj.dimension
                    
                    Dmat = constructMatrix(obj.x{k}, diffDegrees(k, j), accuracy);
                    
                    D{j} = kron(Dmat, D{j});
                end
            end
            
            function mat = constructMatrix(x, degree, accuracy)
                [stencil, coefficients] = centredScheme(degree, accuracy);
                xLength = length(x);
                dx = x(2)-x(1);
                B = repmat(coefficients*dx.^(-degree),xLength,1);
                mat = spdiagsPeriodicBCs(B, stencil, xLength);
            end
            
            function [stencil,coefficients] = centredScheme(degree, accuracy)
                stencilLength = ceil(degree/2) + ceil(accuracy/2)-1;
                stencil = -stencilLength:stencilLength;
                ex = (0:2*stencilLength)';
                S = repmat(stencil,stencilLength*2+1,1).^ex;
                del = factorial(degree) * (ex == degree);
                coefficients = (S\del)';
            end
            
            function mat = spdiagsPeriodicBCs(B, stencil, size1)
                mat = spdiags(B, stencil, size1, size1);
                mat = spdiags(B, stencil + size1, size1, size1) + mat;
                mat = spdiags(B, stencil - size1, size1, size1) + mat;
            end
            
        end
        
        function Y = reshapeToVector(obj, y)
            Y = reshape(y, prod(obj.shape), []);
        end
        
        function index = findDiffCellIndex(obj, deg)
            [~, index] = ismember(deg(:,1)', obj.diffDegrees', 'rows');
        end
        
        function y = reshapeToShape(obj, Y)
            if obj.dimension == 1
                y = reshape(Y, [obj.shape(1), size(Y, 2)]);
            else
                y = reshape(Y, [obj.shape, size(Y, 2)]);
            end
        end
    end
end