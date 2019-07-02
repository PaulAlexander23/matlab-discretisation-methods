function tests = tests()
    tests = functiontests(localfunctions);
end

function testFiniteDifferenceCoefficients(testCase)
    MaxDiff = 4; NumberOfPoints = 8; centrePoint = 0; pointFromCentre = [0,1,-1,2,-2,3,-3,4,-4];
    
    d = finiteDifferenceCoefficients(MaxDiff, NumberOfPoints, centrePoint, pointFromCentre);
    
    expected = [91/8, -122/15, -122/15, 169/60, 169/60, -2/5, -2/5, 7/240, 7/240]';
    
    verifyEqual(testCase, squeeze(d(4 + 2,end,:)), expected, 'AbsTol', 1e-15);
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

function test1dCentralFiniteDifferenceMatrixConstruction(testCase)
    x = {1:16};
    domain = FDDomain(x, 1, 2);
    actual = domain.getDiffMatrix(1);
    B = [-0.5,0.5].*ones(16,1);
    expected = spdiags(B, [-1,1], 16, 16);
    expected(16,1) = 0.5;
    expected(1,16) = -0.5;
    verifyEqual(testCase, actual,expected);
end

function test1dForwardFiniteDifferenceMatrixConstruction(testCase)
    x = {1:16};
    domain = FDDomain(x, 1, 2);
    actual = domain.getDiffMatrix(1);
    B = [-0.5,0.5].*ones(16,1);
    expected = spdiags(B, [-1,1], 16, 16);
    expected(16,1) = 0.5;
    expected(1,16) = -0.5;
    verifyEqual(testCase, actual,expected);
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

function [actual, expected] = function1d(domain)
    y = cos(2*pi*domain.x{1});
    dy = domain.diff(y, 1);
    d2y = domain.diff(y, 2);
    actual = d2y / (4 * pi^2) + dy / (2 * pi) + y;
    expected = -sin(2*pi*domain.x{1});
end

function testEvaluatingFunction1dFiniteDifference(testCase)
    x = setup1dX(2^8);
    
    domain = FDDomain(x, [1, 2], 2);
    
    [actual, expected] = function1d(domain);
    
    verifySize(testCase, actual, [2^8, 1])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
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

function [actual, expected] = function2d(domain)
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');
    
    dy1 = domain.diff(y, [1, 0]');
    dy2 = domain.diff(y, [0, 1]');
    
    actual = dy1 / 2 / pi + dy2 / 2 / pi;
    expected = -sin(2*pi*domain.x{1}) - sin(2*pi*domain.x{2}');
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

function testEvaluatingFunction2dVectorisedFiniteDifference(testCase)
    x = setup2dX(2^8);
    
    problemDeg = [1, 0; 0, 1]';
    
    domain = FDDomain(x, problemDeg, 2);
    
    [actual, expected] = function2dVectorised(domain);
    
    verifySize(testCase, actual, [2^8, 2^8, 2])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test1dBoundaryConditionsFiniteDifference(testCase)
    x = setup1dX(2^8);
    
    c0 = [1, 0];
    c1 = [0, 1];
    f = [0, 0];
    
    domain = FDDomain(x, 1, 4, {[c0, c1, f]});
    
    y0 = domain.x{1} .* (domain.x{1} - 1).^3;
    expected = (domain.x{1} - 1).^3 + 3 * domain.x{1}.^2 .* (domain.x{1} - 1).^2;
    actual = domain.diff(y0, 1);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test1dPeriodicBoundaryConditionsFiniteDifference(testCase)
    x = setup1dX(2^8);
    
    domain = FDDomain(x, 1, 4, {"Periodic"});
    
    y0 = cos(2 * pi * domain.x{1});
    expected = -2 * pi * sin(2 * pi * domain.x{1});
    actual = domain.diff(y0, 1);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test2dBoundaryConditionsFiniteDifference(testCase)
    x = setup2dX(2^8);
    
    xc0 = [ones(2^8,1), zeros(2^8,1)];
    xc1 = [zeros(2^8,1), zeros(2^8,1)];
    xf = [zeros(2^8,1), zeros(2^8,1)];
    yc0 = [ones(2^8,1), zeros(2^8,1)];
    yc1 = [zeros(2^8,1), ones(2^8,1)];
    yf = [zeros(2^8,1), zeros(2^8,1)];
    
    domain = FDDomain(x, [1, 0; 0, 1]', 4, {[xc0, xc1, xf], [yc0, yc1, yf]});
    
    y0 = (domain.x{1}.^2 /2 - domain.x{1}) .* (cos(2*pi * domain.x{2}') - 1);
    expectedx = (cos(2*pi * domain.x{2}') - 1) .* (domain.x{1} - 1);
    expectedy = -2*pi * sin(2*pi * domain.x{2}') .* (domain.x{1}.^2 /2 - domain.x{1});
    
    verifyEqual(testCase, domain.diff(y0, [1; 0]), expectedx, 'RelTol', 1e-3, 'AbsTol', 1e-4)
    verifyEqual(testCase, domain.diff(y0, [0; 1]), expectedy, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test2dPeriodicBoundaryConditionsFiniteDifference1(testCase)
    x = setup2dX(2^8);
    
    xc0 = [ones(2^8,1), zeros(2^8,1)];
    xc1 = [zeros(2^8,1), zeros(2^8,1)];
    xf = [zeros(2^8,1), zeros(2^8,1)];
    
    domain = FDDomain(x, [1, 0; 0, 1]', 4, {[xc0, xc1, xf], "Periodic"});
    
    y0 = (domain.x{1}.^2 /2 - domain.x{1}) .* cos(2*pi * domain.x{2}');
    expectedx = cos(2*pi * domain.x{2}') .* (domain.x{1} - 1);
    expectedy = -2*pi * sin(2*pi * domain.x{2}') .* (domain.x{1}.^2 /2 - domain.x{1});
    
    verifyEqual(testCase, domain.diff(y0, [1; 0]), expectedx, 'RelTol', 1e-3, 'AbsTol', 1e-4)
    verifyEqual(testCase, domain.diff(y0, [0; 1]), expectedy, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test2dPeriodicBoundaryConditionsFiniteDifference2(testCase)
    x = setup2dX(2^8);
    
    domain = FDDomain(x, [1, 0; 0, 1]', 4, {"Periodic", "Periodic"});
    
    y0 = cos(2*pi * domain.x{1}) .* cos(2*pi * domain.x{2}');
    expectedx = cos(2*pi * domain.x{2}') .* -2*pi .* sin(2*pi * domain.x{1});
    expectedy = -2*pi * sin(2*pi * domain.x{2}') .* cos(2*pi * domain.x{1});
    verifyEqual(testCase, domain.diff(y0, [1; 0]), expectedx, 'RelTol', 1e-3, 'AbsTol', 1e-4)
    verifyEqual(testCase, domain.diff(y0, [0; 1]), expectedy, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function x = setup1dX(n)
    x = {linspace(1/n, 1, n)'};
end

function x = setup2dX(n)
    x = {linspace(1/n, 1, n)'; linspace(1/n, 1, n)'};
end
