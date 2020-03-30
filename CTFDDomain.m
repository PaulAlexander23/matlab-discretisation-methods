classdef CTFDDomain < Domain
    properties
        DA
        DB
        diffDegrees
        accuracy
    end
    
    methods
        function obj = CTFDDomain(x, diffDegrees, accuracy)
            obj = obj@Domain(x);
            
            obj.diffDegrees = diffDegrees;
            obj.accuracy = accuracy;
            [obj.DA, obj.DB] = initialiseDiffMatrices(obj, obj.diffDegrees, obj.accuracy);
        end
        
        function dy = diff(obj, y, degree)
            dy = reshapeToVector(obj, y);
            
            [DxA, DxB] = obj.diffMat(degree);
            
            largeDxA = kron(speye(size(dy,1)/size(DxA,1)), DxA);
            largeDxB = kron(speye(size(dy,1)/size(DxB,1)), DxB);
            
            dy = largeDxA \ (largeDxB * dy);
            
            dy = reshapeToDomain(obj, dy);
        end
        
        function [D, D2] = diffMat(obj, degree)
            index = findDiffCellIndex(obj,degree);
            
            if nargout == 2
                D = obj.DA{index};
                D2 = obj.DB{index};
            elseif nargout == 1
                D = obj.DA{index} \ obj.DB{index};
            end
        end
        
        function w = multiply(~, u, v, powers)
            if nargin < 4, powers = [1, 1]; end
            
            w = u.^powers(1) .* v.^powers(2);
        end
    end
    
    methods(Access = private)
        function [DA, DB] = initialiseDiffMatrices(obj, diffDegrees, accuracy)
            DA = cell(1, size(diffDegrees, 2));
            DB = cell(1, size(diffDegrees, 2));
            
            for j = 1:size(diffDegrees, 2)
                DA{j} = 1;
                DB{j} = 1;
                for k = 1:obj.dimension
                    [DAmat, DBmat] = constructMatrices(obj.x{k}, diffDegrees(k, j), accuracy);
                    
                    DA{j} = kron(DAmat, DA{j});
                    DB{j} = kron(DBmat, DB{j});
                end
            end
            
            function [DAmat, DBmat] = constructMatrices(x, degree, accuracy)
                
                lhsStencilLength = 1 + 2 * fix((degree + accuracy - 1)/4);
                rhsStencilLength = 1 + 2 * fix((degree + accuracy + 1)/4);
                
                [lhsStencil, lhsCoefficients, rhsStencil, rhsCoefficients] = ...
                    ctfdScheme(lhsStencilLength, rhsStencilLength, degree);
                
                xLength = length(x);
                dx = x(2) - x(1);
                
                A = repmat(lhsCoefficients.', xLength, 1);
                DAmat = spdiagsPeriodicBCs(A, lhsStencil, xLength);
                
                B = repmat(rhsCoefficients.'*dx.^(-degree), xLength, 1);
                DBmat = spdiagsPeriodicBCs(B, rhsStencil, xLength);
            end
            
            function [lhsStencil, lhsCoefficients, rhsStencil, rhsCoefficients] = ...
                    ctfdScheme(lhsStencilLength, rhsStencilLength, degree, alpha)
                if nargin < 4, alpha = []; end
                
                lhsStencil = (ceil(-(lhsStencilLength-1)/2):ceil((lhsStencilLength-1)/2))';
                rhsStencil = (ceil(-(rhsStencilLength-1)/2):ceil((rhsStencilLength-1)/2))';
                
                totalNumberOfCoefficients = lhsStencilLength + rhsStencilLength;
                
                A = [zeros(lhsStencilLength, degree), ...
                    taylorSeries(lhsStencil', totalNumberOfCoefficients - 1 - degree)'];
                B = taylorSeries(rhsStencil', totalNumberOfCoefficients - 1)';
                
                
                M = [A; B].';
                
                lhsCondition = full(ind2vec(find(lhsStencil == 0), totalNumberOfCoefficients)).';
                
                M = [M; lhsCondition];
                
                r = full(ind2vec(totalNumberOfCoefficients));
                
                for m = 1:length(alpha)
                    ebm = ind2vec(find(lhsStencil == m), totalNumberOfCoefficients).';
                    ebminusm = ind2vec(find(lhsStencil == -m), totalNumberOfCoefficients).';
                    
                    M(end-m, :) = ebminusm.' - alpha(m) * ebm.';
                end
                
                x = M\r;
                
                lhsCoefficients = x(1:lhsStencilLength);
                rhsCoefficients = -x(lhsStencilLength+1:lhsStencilLength + rhsStencilLength);
                
                function series = taylorSeries(dx, numberOfTerms)
                    n = (0:numberOfTerms-1)';
                    series = dx .^ n ./factorial(n);
                end
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
