function F = fbenney(x,y,params,method)
    deg = [1,0;2,0;0,2]';
    delta = params(1);
    theta = params(2);
    Re = params(3);
    We = params(4);
    C = params(5);
    
    L = cellfun(@(x) x(end),x);
    e1 = zeros(1,1,2);
    e1(1,1,1) = 1;
    dy = feval(method,x,y,deg);
    
    P = y * cot(theta) - 1/2 * 1/C * (dy(:,2)+dy(:,3));
    %P = H * cot(theta) - We * R(H, L) - 1/2 * 1/C * lap(H, L);
    
    F = div(x,...
        (2/3 * y.^3 + 8/15 * Re * delta * y .^ 6 .* dy(:,1)) .* e1...
        - 2/3 * delta * y.^3 .* grad(x,P));
    
    function out = div(x,y)
        divdeg = cat(3,[1,0;0,0],[0,0;0,1]);
        d2y = diff_ps_2d(x,y,divdeg);
        out = d2y(:,:,1) + d2y(:,:,2);
    end
    
    function out = grad(x,y)
        graddeg = [1,0;0,1]';
        out = diff_ps_2d(x,y,graddeg);
    end
end