classdef PSDomain < Domain
    properties
        length
        wavenumber
        suppression
        normaliseAmplitude
    end
    
    methods
        function obj = PSDomain(x)
            obj = obj@Domain(x);
            obj.length = calculateLength(obj);
            obj.wavenumber = calculateWavenumber(obj);
            obj.suppression = 1e-13;
            obj.normaliseAmplitude = (2/prod(obj.shape)) / obj.fftScaling();
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
            
            upad = obj.ifftn(obj.zeropad(uhat,ratio)) * ratio.^obj.dimension;
            vpad = obj.ifftn(obj.zeropad(vhat,ratio)) * ratio.^obj.dimension;
            
            wpad = upad.^powers(1) .* vpad.^powers(2);
            
            what = obj.trunc( ...
                obj.fftn(wpad) / ratio.^obj.dimension ...
                , 1/ratio);
        end
        
        function f = fftn(obj, x)
            f = fftn(x) * obj.fftScaling();
        end
        
        function f = ifftn(obj, x)
            f = ifftn(x, 'symmetric') / obj.fftScaling();
        end

        function uhatpad = zeropad(obj, uhat, ratio)
            if obj.dimension == 1
                uhatpad = obj.zeropad1d(uhat,ratio);
            elseif obj.dimension == 2
                uhatpad = permute(obj.zeropad1d(permute(obj.zeropad1d(uhat,ratio),[2,1,3]),ratio),[2,1,3]);
            else
                error("zeropad not defined for higher dimensions.")
            end
        end
        
        function uhat = trunc(obj, uhatpad, ratio)
            if obj.dimension == 1
                uhat = obj.trunc1d(uhatpad,ratio);
            elseif obj.dimension == 2
                uhat = permute(obj.trunc1d(permute(obj.trunc1d(uhatpad,ratio),[2,1,3]),ratio),[2,1,3]);
            else
                error("trunc not defined for higher dimensions.")
            end
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
        
        function scaling = fftScaling(obj)
            if obj.dimension == 1
                scaling = obj.length / obj.shape(1);
            else
                scaling = prod(obj.length ./ obj.shape');
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
        
        function upad = zeropad1d(obj, u, ratio)
            s = size(u);
            upad = zeros([s(1) * ratio, s(2:end)]);
            ind = [1:s(1)/2, s(1) * (ratio - 0.5) + 1:s(1) * ratio];
            upad(ind,:,:) = u;
        end
        
        function u = trunc1d(obj, upad, ratio)
            N = size(upad,1);
            ind = [1:ratio * N/2, N * (1 - ratio/2) + 1:N];
            u = upad(ind,:,:);
        end
    end
end