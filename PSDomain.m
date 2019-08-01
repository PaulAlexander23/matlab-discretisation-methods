classdef PSDomain < Domain
    properties
        length
        wavenumber
        suppression
    end
    
    methods
        function obj = PSDomain(x)
            obj = obj@Domain(x);
            obj.length = calculateLength(obj);
            obj.wavenumber = calculateWavenumber(obj);
            obj.suppression = 1e-13;
        end
        
        function dyhat = diff(obj, yhat, degree)
            dyhat = priorSuppression(obj, yhat);
            
            dyhat = wavenumberMultiplicand(obj, degree).*dyhat;
        end

        function what = multiply(obj, uhat, vhat, powers) % Convolution Really
            if nargin < 4
                powers = [1, 1];
            end
            ratio = (sum(abs(powers)) + 1)/2;
            
            upad = obj.ifftn(obj.matrixZeropad(uhat,ratio)) * ratio.^obj.dimension;
            vpad = obj.ifftn(obj.matrixZeropad(vhat,ratio)) * ratio.^obj.dimension;
            
            wpad = upad.^powers(1) .* vpad.^powers(2);
            
            what = obj.matrixTrunc( ...
                obj.fftn(wpad) / ratio.^obj.dimension ...
                , ratio);
        end
        
        function f = fftn(obj, x)
            dt = prod(obj.length ./ obj.shape');
            f = fftn(x) * dt;
        end
        
        function f = ifftn(obj, x)
            dt = prod(obj.length ./ obj.shape');
            f = ifftn(x, 'symmetric') / dt;
        end
    end
    
    methods (Access=private)
        function L = calculateLength(obj)
            if obj.x{1}(1) ~= 0
                L = cellfun(@(x) x(end), obj.x);
            else
                L = cellfun(@(x) x(end) + x(1), obj.x);
            end
        end
        
        function k = calculateWavenumber(obj)
            k = cell(1, obj.dimension);
            for d = 1:obj.dimension
                k{d} = [0:obj.shape(d)/2 - 1, 0, 1 - obj.shape(d)/2:-1]' * ...
                    2*pi/obj.length(d);
            end
        end
        
        function dyhat = priorSuppression(obj, dyhat)
            dyhat(abs(dyhat)<obj.suppression) = 0;
        end
        
        function f = wavenumberMultiplicand(obj, degree)
            if obj.dimension == 1
                f = (1i*obj.wavenumber{1}).^degree;
            elseif obj.dimension == 2
                f = (1i*obj.wavenumber{1}).^degree(1)*(1i*obj.wavenumber{2}').^degree(2);
            else
                error("diff not defined for higher dimensions.")
            end
        end
        
        function upad = zeropad(obj, u, ratio)
            s = size(u);
            upad = zeros([s(1) * ratio, s(2:end)]);
            ind = [1:s(1)/2, s(1) * (ratio - 0.5) + 1:s(1) * ratio];
            upad(ind,:,:) = u;
        end
        
        function uhatpad = matrixZeropad(obj, uhat, ratio)
            if obj.dimension == 1
                uhatpad = obj.zeropad(uhat,ratio);
            elseif obj.dimension == 2
                uhatpad = permute(obj.zeropad(permute(obj.zeropad(uhat,ratio),[2,1,3]),ratio),[2,1,3]);
            else
                error("matrixZeropad not defined for higher dimensions.")
            end
        end
        
        function u = trunc(obj, upad, ratio)
            N = size(upad,1);
            ind = [1:ratio * N/2, N * (1 - ratio/2) + 1:N];
            u = upad(ind,:,:);
        end
        
        function uhat = matrixTrunc(obj, uhatpad, ratio)
            if obj.dimension == 1
                uhat = obj.trunc(uhatpad,1/ratio);
            elseif obj.dimension == 2
                uhat = permute(obj.trunc(permute(obj.trunc(uhatpad,1/ratio),[2,1,3]),1/ratio),[2,1,3]);
            else
                error("matrixTrunc not defined for higher dimensions.")
            end
        end
    end
end