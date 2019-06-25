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

        function dy = diff(obj, y, degree)
            y = reshapeToDomain(obj, y);

            dy = fftn(y);

            dy = priorSuppression(obj, dy);

            dy = wavenumberMultiplicand(obj, degree) .* dy;

            dy = ifftn(dy);

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
                k{d} = [0:obj.shape(d) / 2 - 1, 0, 1 - obj.shape(d) / 2:-1]' * ...
                    2 * pi / obj.length(d);
            end
        end

        function dy = priorSuppression(obj, dy)
            dy(abs(dy) < obj.suppression) = 0;
        end

        function f = wavenumberMultiplicand(obj, degree)
            if obj.dimension == 1
                f = (1i * obj.wavenumber{1}).^degree;
            elseif obj.dimension == 2
                f = (1i * obj.wavenumber{1}).^degree(1) * (1i * obj.wavenumber{2}').^degree(2);
            else
                error("diff not defined for higher dimensions.")
            end
        end
    end
end