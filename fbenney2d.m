function F = fbenney2d(domain, y, params)
    
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);
    
    e1 = zeros(1, 1, 2);
    e1(1, 1, 1) = 1;
    
    P = y * cot(theta) - 1/2 * 1 / C * (domain.diff(y, [2; 0]) + domain.diff(y, [0; 2]));
    
    gradP = cat(3, domain.diff(P, [1; 0]), domain.diff(P, [0; 1]));
    
    Q = (2/3 * domain.multiply(y, y, [2, 1]) + 8/15 * Re * delta * domain.multiply(y, domain.diff(y, [1; 0]), [6, 1])) .* e1 ...
        -2/3 * delta * domain.multiply(y, gradP, [3, 1]);
    
    F = - domain.diff(Q(:, :, 1), [1, 0]') - ...
        domain.diff(Q(:, :, 2), [0, 1]');
    
end
