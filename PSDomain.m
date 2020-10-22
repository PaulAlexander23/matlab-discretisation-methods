classdef PSDomain < Domain
    properties
        length
        wavenumber
        suppression
        antialiasing
        scaling
        complex
    end

    properties (Access = private)
        fourierDomain
    end

    methods
        function obj = PSDomain(x, antialiasing, complex, suppression)
            if nargin < 2, antialiasing = false; end
            if nargin < 3, complex = true; end
            if nargin < 4, suppression = 1e-15; end
            
            obj = obj@Domain(x);
            obj.length = calculateLength(obj);
            obj.suppression = suppression;
            obj.antialiasing = antialiasing;
            obj.scaling = 2/prod(obj.shape);
            obj.complex = complex;
            obj.wavenumber = calculateWavenumber(obj);
            obj.fourierDomain = Domain(obj.wavenumber);
        end
        
        function Y = reshapeToVector(obj, y)
            Y = obj.fourierDomain.reshapeToVector(y);
        end
        
        function y = reshapeToDomain(obj, Y)
            y = obj.fourierDomain.reshapeToDomain(Y);
        end
        
        function f = fft(obj, x)
            xcell = obj.extractSurfacesAsCell(x);
            fcell = obj.fftCell(xcell);
            f = cell2mat(fcell);
            
            f = zeroSmallModes(obj, f);
        end
        
        function x = ifft(obj, f)
            fcell = obj.fourierDomain.extractSurfacesAsCell(f);
            xcell = obj.ifftCell(fcell);
            x = cell2mat(xcell);
        end
        
        function dyhat = diff(obj, yhat, degree)
            yhatCell = obj.fourierDomain.extractSurfacesAsCell(yhat);
            dyhatCell = obj.diffCell(yhatCell, degree);
            dyhat = cell2mat(dyhatCell);
        end
        
        function D = diffMat(obj, degree)
            D = wavenumberMultiplicand(obj, degree);
        end
        
        function what = multiply(obj, uhat, vhat, powers)
            if nargin < 4, powers = [1, 1]; end
            
            uhatCell = obj.fourierDomain.extractSurfacesAsCell(uhat);
            vhatCell = obj.fourierDomain.extractSurfacesAsCell(vhat);
            whatCell = obj.multiplyCell(uhatCell, vhatCell, powers);
            what = cell2mat(whatCell);
        end

        function uhat = zeroSmallModes(obj, uhat, tolerance)
            if nargin < 3, tolerance = obj.suppression; end

           uhat(abs(uhat) < tolerance) = 0; 
        end

        function out = filterOutShortWaves(obj, in, ratioKeptToAll, mask)
            if nargin < 4, mask = "rectangle"; end

            ratio = ratioKeptToAll/2;
            if obj.dimension == 1
                shortWaveFilter = ones(obj.fourierDomain.shape);
                shortWaveFilter(floor(end*ratio+1):ceil(end-end*ratio)) = 0;
            elseif obj.dimension == 2
                shortWaveFilter = zeros(obj.fourierDomain.shape);
                for k = 1:obj.fourierDomain.shape(2)
                    for j = 1:obj.fourierDomain.shape(1)
                        if inRegion(j,k,obj.shape,ratio .* obj.shape, mask)
                            shortWaveFilter(j,k) = 1;
                        end
                    end
                end
            end
            
            out = in .* shortWaveFilter;
            
            function out = inRegion(j,k,shape,radii, mask)
                if mask == "ellipse"
                    out = any([ ...
                        inEllipse(j,k,[1,1],radii), ...
                        inEllipse(j,k,[shape(1),1],radii), ...
                        inEllipse(j,k,[1,shape(2)],radii), ...
                        inEllipse(j,k,[shape(1),shape(2)],radii)]);
                elseif mask == "triangle"
                    out = any([ ...
                        inTriangle(j,k,[1,1],radii), ...
                        inTriangle(j,k,[shape(1),1],radii), ...
                        inTriangle(j,k,[1,shape(2)],radii), ...
                        inTriangle(j,k,[shape(1),shape(2)],radii)]);
                elseif mask == "rectangle"
                    out = any([ ...
                        inRectangle(j,k,[1,1],radii), ...
                        inRectangle(j,k,[shape(1),1],radii), ...
                        inRectangle(j,k,[1,shape(2)],radii), ...
                        inRectangle(j,k,[shape(1),shape(2)],radii)]);
                    
                end
            end
            
            function out = inEllipse(j,k,centre,radii)
                out = (j - centre(1)).^2 ./ radii(1).^2 + (k - centre(2)).^2 ./ radii(2).^2 < 1;
                
            end
            
            function out = inTriangle(j, k, centre, sides)
                out = abs(j - centre(1)) ./ sides(1) + abs(k - centre(2)) ./ sides(2) < 1;
            end

            function out = inRectangle(j, k, centre, sides)
                out = (abs(j - centre(1)) ./ sides(1)) < 1 && (abs(k - centre(2)) ./ sides(2)) < 1;
            end
        end
    end

    methods(Access = private)
        function L = calculateLength(obj)
            L = cellfun(@(x) 2*x(end) - x(end-1) - x(1), obj.x);
        end

        function k = calculateWavenumber(obj)
            k = cell(1, obj.dimension);
            for d = 1:obj.dimension
                if ~obj.complex && d == 1
                    w = (0:obj.shape(d)/2 - 1)';
                else
                    w = [0:obj.shape(d)/2 - 1, 0, 1 - obj.shape(d)/2:-1]';
                end
                
                k{d} = w * (2*pi/obj.length(d));
            end
        end

        function f = removeSymmetricConjugate(~, f)
            f = f(1:end/2,:,:);
        end
        
        function f = addSymmetricConjugate(~, f)
            f = [f; zeros(size(f))];
        end

        function fcell = fftCell(obj, xcell)
            fcell = cellfun(@(x)obj.fftLocal(x), xcell, 'UniformOutput', false);
        end

        function f = fftLocal(obj, x)
            f = fftn(x) * obj.scaling;

            if ~obj.complex
                f = obj.removeSymmetricConjugate(f);
            end
        end

        function xcell = ifftCell(obj, fcell)
            xcell = cellfun(@(f) obj.ifftLocal(f), fcell, 'UniformOutput', false);
        end

        function x = ifftLocal(obj, f)
            if obj.complex
                x = ifftn(f) / obj.scaling;
            else
                f = obj.addSymmetricConjugate(f);
                x = ifftn(f, 'symmetric') / obj.scaling;
            end
            
        end

        function dycell = diffCell(obj, ycell, degree)
            dycell = cellfun(@(y) obj.diffLocal(y, degree), ycell, 'UniformOutput', false);
        end

        function dyhat = diffLocal(obj, yhat, degree)
            dyhat = obj.zeroSmallModes(yhat);
            
            dyhat = obj.reshapeToVector(dyhat);
            
            dyhat = obj.diffMat(degree) * dyhat;
            
            dyhat = obj.reshapeToDomain(dyhat);
        end

        function whatCell = multiplyCell(obj, uhatCell, vhatCell, powers)
            [M, N] = size(uhatCell);
            whatCell = cell(M, N);
            for m = 1:M
                for n = 1:N
                    whatCell{m,n} = obj.multiplyLocal(uhatCell{m,n}, vhatCell{m,n}, powers);
                end
            end
        end

        function what = multiplyLocal(obj, uhat, vhat, powers)
            if obj.antialiasing
                ratio = (sum(abs(powers)) + 1)/2;

                uhat = obj.zeropad(uhat,ratio) * ratio.^obj.dimension;
                vhat = obj.zeropad(vhat,ratio) * ratio.^obj.dimension;
            end

            u = obj.ifftLocal(uhat);
            v = obj.ifftLocal(vhat);

            w = u.^powers(1) .* v.^powers(2);
            
            what = obj.fftLocal(w);
            
            if obj.antialiasing
                what = obj.trunc(what / ratio.^obj.dimension, ...
                    1/ratio);
            end
        end

        function uhatpad = zeropad(obj, uhat, ratio)
            if obj.dimension == 1
                uhatpad = obj.zeropad1d(uhat, ratio);
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

        function f = wavenumberMultiplicand(obj, degree)
            if obj.complex
                m = prod(obj.shape);
            else
                m = prod(obj.shape)/2;
            end
            if obj.dimension == 1
                A = (1i * obj.wavenumber{1}).^degree;
            elseif obj.dimension == 2
                A = ((1i*obj.wavenumber{1}).^degree(1) * ...
                    (1i*obj.wavenumber{2}').^degree(2));
                A = reshapeToVector(obj, A);
            else
                error("diff not defined for higher dimensions.")
            end
            f = spdiags(A,0,m,m);
        end
        
        function upad = zeropad1d(obj, u, ratio, dimension)
            if nargin < 4, dimension = 1; end
            
            s = size(u);
            N = s(1) * ratio;
            upad = zeros([N, s(2:end)]);
            if ~obj.complex && dimension == 1
                ind = 1:s(1);
            else
                ind = [1:s(1)/2, N - s(1)/2 + 1:N];
            end
            upad(ind,:,:) = u;
        end
        
        function u = trunc1d(obj, upad, ratio, dimension)
            if nargin < 4, dimension = 1; end
            
            N = size(upad,1);
            s = size(upad);
            u = zeros([N*ratio, s(2:end)]); 
            if ~obj.complex && dimension == 1
                ind = (1:ratio * N)';
                u(ind,:,:) = upad(ind,:,:);
            else
                ind = [1:ratio * N/2, N - ratio * N/2 + 2:N]';
                u([1:ratio*N/2,ratio*N/2+2:ratio*N],:,:) = upad(ind,:,:);
            end
        end
    end
end
