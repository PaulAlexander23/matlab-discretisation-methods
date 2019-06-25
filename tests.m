function tests = tests()
    tests = functiontests(localfunctions);
end

function testDomain1d(testCase)
    expectedX = setup1dX(2^8);
    expectedShape = [2^8, 1];

    domain = Domain(expectedX);
    actualX = domain.x;
    actualShape = domain.shape;

    verifyEqual(testCase, actualX, expectedX)
    verifyEqual(testCase, actualShape, expectedShape)
    verifyEqual(testCase, domain.dimension, 1)
end

function testDomain2d(testCase)
    expectedX = setup2dX(2^8);
    expectedShape = [2^8, 2^8];

    domain = Domain(expectedX);
    actualX = domain.x;
    actualShape = domain.shape;

    verifyEqual(testCase, actualX, expectedX)
    verifyEqual(testCase, actualShape, expectedShape)
    verifyEqual(testCase, domain.dimension, 2)
end

function test1dFiniteDifference(testCase)
    x = setup1dX(2^8);

    domain = FDDomain(x, 1, 2);

    Y = cos(2*pi*domain.x{1});

    degree = 1;

    actual = domain.diff(Y, degree);
    expected = -2 * pi * sin(2*pi*x{1});

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test1dPseudoSpectral(testCase)
    x = setup1dX(2^8);

    domain = PSDomain(x);

    verifyEqual(testCase, domain.length, 1);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        [0:2^7 - 1, 0, 1 - 2^7:-1]'*2*pi);

    Y = cos(2*pi*domain.x{1});

    degree = 1;

    actual = domain.diff(Y, degree);
    expected = -2 * pi * sin(2*pi*x{1});

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14)
end

function test1dFiniteDifferenceGetDiffMatrix(testCase)
    x = setup1dX(2^8);

    domain = FDDomain(x, 1, 2);

    degree = 1;

    actual = domain.getDiffMatrix(degree);
    expectedSize = [2^8, 2^8];

    verifySize(testCase, actual, expectedSize);
end

function testEvaluatingFunction1dFiniteDifference(testCase)
    x = setup1dX(2^8);

    domain = FDDomain(x, [1, 2], 2);

    [actual, expected] = function1d(domain);

    verifySize(testCase, actual, [2^8, 1])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function testEvaluatingFunction1dVectorisedFiniteDifference(testCase)
    x = setup1dX(2^8);

    domain = FDDomain(x, [1, 2], 2);

    [actual, expected] = function1dVectorised(domain);

    verifySize(testCase, actual, [2^8, 2])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function testEvaluatingFunction1dPseudoSpectral(testCase)
    x = setup1dX(2^8);

    domain = PSDomain(x);

    [actual, expected] = function1d(domain);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-15)
end

function test2dFiniteDifference(testCase)
    x = setup2dX(2^8);

    domain = FDDomain(x, [1, 0]', 2);

    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');

    degree = [1, 0]';

    actual = domain.diff(Y, degree);
    expected = -2 * pi * sin(2*pi*x{1}) .* ones(size(x{2}'));

    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test2dPseudoSpectral(testCase)
    x = setup2dX(2^8);

    domain = PSDomain(x);

    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');

    degree = [1, 0]';

    actual = domain.diff(Y, degree);
    expected = -2 * pi * sin(2*pi*x{1}) .* ones(size(x{2}'));

    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function testEvaluatingFunction2dFiniteDifference(testCase)
    x = setup2dX(2^8);

    problemDeg = [1, 0; 0, 1]';

    domain = FDDomain(x, problemDeg, 2);

    [actual, expected] = function2d(domain);

    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function testEvaluatingFunction2dPseudoSpectral(testCase)
    x = setup2dX(2^8);

    domain = PSDomain(x);

    [actual, expected] = function2d(domain);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-13)
end

function testEvaluatingFunction2dVectorisedFiniteDifference(testCase)
    x = setup2dX(2^8);

    problemDeg = [1, 0; 0, 1]';

    domain = FDDomain(x, problemDeg, 2);

    [actual, expected] = function2dVectorised(domain);

    verifySize(testCase, actual, [2^8, 2^8, 2])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function [actual, expected] = function1d(domain)
    y = cos(2*pi*domain.x{1});
    dy = domain.diff(y, 1);
    d2y = domain.diff(y, 2);
    actual = d2y / (4 * pi^2) + dy / (2 * pi) + y;
    expected = -sin(2*pi*domain.x{1});
end

function [actual, expected] = function1dVectorised(domain)
    y = cos(2*pi*domain.x{1});
    y2 = cos(2*pi*domain.x{1});
    dy = domain.diff([y, y2], 1);
    d2y = domain.diff([y, y2], 2);
    actual = d2y / (4 * pi^2) + dy / (2 * pi) + y;
    expected = cat(2, ...
        -sin(2*pi*domain.x{1}), ...
        -sin(2*pi*domain.x{1}));
end

function [actual, expected] = function2d(domain)
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');

    dy1 = domain.diff(y, [1, 0]');
    dy2 = domain.diff(y, [0, 1]');

    actual = dy1 / 2 / pi + dy2 / 2 / pi;
    expected = -sin(2*pi*domain.x{1}) - sin(2*pi*domain.x{2}');
end

function [actual, expected] = function2dVectorised(domain)
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');
    y2 = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');

    dy1 = domain.diff([y, y2], [1, 0]');
    dy2 = domain.diff([y, y2], [0, 1]');

    actual = dy1 / 2 / pi + dy2 / 2 / pi;
    expected = cat(3, ...
        -sin(2*pi*domain.x{1})-sin(2*pi*domain.x{2}'), ...
        -sin(2*pi*domain.x{1})-sin(2*pi*domain.x{2}'));
end

function x = setup1dX(n)
    x = {linspace(1/n, 1, n)'};
end

function x = setup2dX(n)
    x = {linspace(1/n, 1, n)'; linspace(1/n, 1, n)'};
end
