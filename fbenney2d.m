function F = fbenney2d(domain, y, params)
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);

    e1 = zeros(1, 1, 2);
    e1(1, 1, 1) = 1;

    P = y * cot(theta) - (domain.diff(y, [2; 0]) + domain.diff(y, [0; 2])) / (2 * C);

    gradP = cat(3, domain.diff(P, [1; 0]), domain.diff(P, [0; 1]));

    Q = (2 / 3 * y.^3 + 8 / 15 * Re * delta * y.^6 .* domain.diff(y, [1; 0])) .* e1 ...
        -2 / 3 * delta * y.^3 .* gradP;

    F = -domain.diff(Q(:, :, 1), [1, 0]') - ...
        domain.diff(Q(:, :, 2), [0, 1]');
end