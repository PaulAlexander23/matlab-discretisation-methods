function dYdt = fwibl1(domain, Y, params)
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);

    y = Y(1:end/2, :, :);
    F1 = Y(end/2+1:end, :, :);

    P = 2 * y * cot(theta) - 1 / C * (domain.diff(y, [2; 0]) + domain.diff(y, [0; 2]));

    F2 = -delta * y.^3 .* domain.diff(P, [0; 1]) / 3;

    dydt = -domain.diff(F1, [1; 0]) - domain.diff(F2, [0; 1]);

    dF1dt = (9 * domain.diff(y, [1; 0]) .* F1.^2 ./ y - ...
        17 * F1 .* domain.diff(F1, [1; 0])) ./ (7 * y) + ...
        (10 * y - 15 * F1 ./ y.^2 - ...
        5 * y .* domain.diff(P, [1; 0])) / (6 * Re);

    dYdt = cat(1, dydt, dF1dt);
end