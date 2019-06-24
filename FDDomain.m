classdef FDDomain < Domain
    properties
        D
        diffDegrees
    end
    
    methods
        function obj = FDDomain(x, diffDegrees, accuracy)
            obj = obj@Domain(x);
            obj.D = init_fd(x, diffDegrees, accuracy);
            obj.diffDegrees = diffDegrees;
            
            initialiseDiffMatrices(obj, diffDegrees, accuracy);
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
        function Y = reshapeToVector(obj, y)
            Y = reshape(y,prod(obj.shape),[]);
        end
        
        function index = findDiffCellIndex(obj, deg)
            [~,index] = ismember(deg(:,1)',obj.diffDegrees','rows');
        end
        
        function y = reshapeToShape(obj, Y)
            y = reshape(Y, [obj.shape,size(Y,2)]);
        end
        
        function D = initialiseDiffMatrices(obj, diffDegrees, accuracy)
            
            dx = cellfun(@(x) x(2)-x(1),obj.x);
            
            D = cell(1,size(diffDegrees,2));
            
            for j = 1:size(diffDegrees,2)
                D{j} = 1;
                for k = 1:obj.dimension
                    [s,c] = centred_scheme(diffDegrees(k,j), accuracy);
                    
                    Dmat = construct_matrix(obj.shape(k),dx(k),s,c,diffDegrees(k,j));
                    
                    D{j} = kron(Dmat,D{j});
                end
            end
            
            function [s,c] = centred_scheme(degree, accuracy)
                sL = ceil(degree/2) + ceil(accuracy/2)-1;
                s = -sL:sL;
                ex = (0:2*sL)';
                S = repmat(s,sL*2+1,1).^ex;
                del = factorial(degree) * (ex == degree);
                c = (S\del)';
            end
            
            function mat = construct_matrix(nx,dx,s,c,d)
                mat = spdiags(ones(nx,1)*c,s,nx,nx)*dx.^(-d);
                mat = mat + spdiags(ones(nx,1)*c,s+nx,nx,nx)*dx.^(-d);
                mat = mat + spdiags(ones(nx,1)*c,s-nx,nx,nx)*dx.^(-d);
            end
        end
    end
end