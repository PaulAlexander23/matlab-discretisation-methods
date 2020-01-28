classdef PSDomain < Domain
    properties
        length
        wavenumber
        suppression
        antialiasing
        scaling
        complex
    end

    methods
        function obj = PSDomain(x, antialiasing, complex, suppression)
            if nargin < 2, antialiasing = false; end
            if nargin < 3, complex = true; end
            if nargin < 4, suppression = eps; end
            
            obj = obj@Domain(x);
            obj.length = calculateLength(obj);
            obj.suppression = suppression;
            obj.antialiasing = antialiasing;
            obj.scaling = 2/prod(obj.shape);
            obj.complex = complex;
            obj.wavenumber = calculateWavenumber(obj);
        end
        
        function Y = reshapeToVector(obj, y)
            if obj.complex
                Y = reshapeToVector@Domain(obj, y);
            else
                Y = reshape(y, prod(obj.shape)/2, []);
            end
        end
        
        function y = reshapeToDomain(obj, Y)
            if obj.complex
                y = reshapeToDomain@Domain(obj, Y);
            else
                if obj.dimension == 1
                    y = reshape(Y, [obj.shape(1)/2, numel(Y) / obj.shape(1)]);
                else
                    y = reshape(Y, [obj.shape(1)/2, obj.shape(2:end), numel(Y) / prod(obj.shape) * 2]);
                end
            end
        end
        
        function dyhat = diff(obj, yhat, degree)
            dyhat = suppress(obj, yhat);
            
            dyhat = reshapeToVector(obj, dyhat);
            
            dyhat = obj.diffMat(degree) * dyhat;
            
            dyhat = reshapeToDomain(obj, dyhat);
        end
        
        function D = diffMat(obj, degree)
            D = wavenumberMultiplicand(obj, degree);
        end
        
        function what = multiply(obj, uhat, vhat, powers)
            if nargin < 4, powers = [1, 1]; end
            
            if obj.antialiasing
                ratio = (sum(abs(powers)) + 1)/2;
                
                uhat = obj.zeropad(uhat,ratio) * ratio.^obj.dimension;
                vhat = obj.zeropad(vhat,ratio) * ratio.^obj.dimension;
            end
            
            u = obj.ifft(uhat);
            v = obj.ifft(vhat);
            
            w = u.^powers(1) .* v.^powers(2);
            
            what = obj.fft(w);
            
            if obj.antialiasing
                what = obj.trunc(what / ratio.^obj.dimension, ...
                    1/ratio);
            end
        end
        
        function f = fft(obj, x)
            if obj.dimension == 1
                f = fft(x);
            elseif obj.dimension == 2
                f = fft2(x);
            else
                error('fft in 3 dimensions and higher is not defined yet.')
            end
            
            if ~obj.complex
                f = f(1:end/2,:,:);
            end
            
            f = f * obj.scaling;
            f = suppress(obj, f);
        end
        
        function x = ifft(obj, f)
            if ~obj.complex
                f = [f; zeros(size(f))];
            end
            
            if obj.dimension == 1
                x = ifft(f, 'symmetric');
            elseif obj.dimension == 2
                x = ifft2(f, 'symmetric');
            else
                error('ifft in 3 dimensions and higher is not defined yet.')
            end
            x = x / obj.scaling;
        end
        
        function uhatpad = zeropad(obj, uhat, ratio)
            if obj.dimension == 1
                uhatpad = obj.zeropad1d(uhat,ratio);
            elseif obj.dimension == 2
                uhatpad = permute(obj.zeropad1d(permute(obj.zeropad1d(uhat,ratio),[2,1,3]),ratio,2),[2,1,3]);
            else
                error("zeropad not defined for higher dimensions.")
            end
        end
        
        function uhat = trunc(obj, uhatpad, ratio)
            if obj.dimension == 1
                uhat = obj.trunc1d(uhatpad,ratio);
            elseif obj.dimension == 2
                uhat = permute(obj.trunc1d(permute(obj.trunc1d(uhatpad,ratio),[2,1,3]),ratio,2),[2,1,3]);
            else
                error("trunc not defined for higher dimensions.")
            end
        end
    end

    methods(Access = private)
        function L = calculateLength(obj)
            if obj.x{1}(1) ~= 0
                L = cellfun(@(x) x(end), obj.x);
            else
                L = cellfun(@(x) x(end)+x(1), obj.x);
            end
        end

        function k = calculateWavenumber(obj)
            k = cell(1, obj.dimension);
            for d = 1:obj.dimension
                if ~obj.complex && d == 1
                    w = (0:obj.shape(d)/2 - 1)';
                else
                    w = [0:obj.shape(d)/2 - 1, 0, 1 - obj.shape(d)/2:-1]';
                end
                
                k{d} = w * ...
                    2*pi/obj.length(d);
            end
        end
        
        function dyhat = suppress(obj, dyhat)
            dyhat(abs(dyhat)<obj.suppression * max(max(abs(dyhat)))) = 0;
        end
        
        function out = filterOutShortWaves(obj, in, ratioKeptToAll)
            ratio = ratioKeptToAll/2;
            if obj.dimension == 1
                f = ones(obj.shape);
                f(round(end*ratio+1):round(end-end*ratio)) = 0;
            elseif obj.dimension == 2
                f = zeros(obj.shape);
                for k = 1:obj.shape(2)
                    for j = 1:obj.shape(1)
                        if inRegion(j,k,obj.shape,ratio .* obj.shape)
                            f(j,k) = 1;
                        end
                    end
                end
            end
            
            out = in .* f;
            
            function out = inRegion(j,k,shape,radii)
                out = any([ ...
                    inEllipse(j,k,[1,1],radii), ...
                    inEllipse(j,k,[shape(1),1],radii), ...
                    inEllipse(j,k,[1,shape(2)],radii), ...
                    inEllipse(j,k,[shape(1),shape(2)],radii)]);
            end
            
            function out = inEllipse(j,k,centre,radii)
                out = (j - centre(1)).^2 ./ radii(1).^2 + (k - centre(2)).^2 ./ radii(2).^2 < 1;
                
            end
            
            function out = inTriangle(j, k, centre, sides)
                out = abs(j - centre(1)) ./ sides(1) + abs(k - centre(2)) ./ sides(2) < 1;
            end
        end
        
        function f = wavenumberMultiplicand(obj, degree)
            if obj.dimension == 1
                A = (1i * obj.wavenumber{1}).^degree;
                f = spdiags(A,0,prod(obj.shape),prod(obj.shape));
            elseif obj.dimension == 2
                A = ((1i*obj.wavenumber{1}).^degree(1) * ...
                    (1i*obj.wavenumber{2}').^degree(2));
                A = reshapeToVector(obj, A);
                
                if obj.complex
                    m = prod(obj.shape);
                else
                    m = prod(obj.shape)/2;
                end
                
                f = spdiags(A,0,m,m);
            else
                error("diff not defined for higher dimensions.")
            end
        end
        
        function upad = zeropad1d(obj, u, ratio, dimension)
            if nargin < 4, dimension = 1; end
            
            s = size(u);
            upad = zeros([s(1) * ratio, s(2:end)]);
            if ~obj.complex && dimension == 1
                ind = 1:s(1);
            else
                ind = [1:s(1)/2, s(1) * (ratio - 0.5) + 1:s(1) * ratio];
            end
            upad(ind,:,:) = u;
        end
        
        function u = trunc1d(obj, upad, ratio, dimension)
            if nargin < 4, dimension = 1; end
            
            N = size(upad,1);
            if ~obj.complex && dimension == 1
                ind = 1:ratio * N;
            else
                ind = [1:ratio * N/2, N * (1 - ratio/2) + 1:N];
            end
            u = upad(ind,:,:);
        end
    end
end