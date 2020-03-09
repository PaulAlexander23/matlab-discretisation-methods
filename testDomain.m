function tests = testDomain()
    tests = functiontests(localfunctions);
end

%% Domain
function testDomain1dConstructor(testCase)
    n = 256;
    x = {linspace(1/n, 1, n)'};
    domain = Domain(x);

    verifyEqual(testCase, domain.x, {linspace(1/256,1,256)'})
    verifyEqual(testCase, domain.shape, [2^8, 1])
    verifyEqual(testCase, domain.dimension, 1)
end

function testDomain2dConstructor(testCase)
    m = 256;
    n = 256;
    x = {linspace(1/m,1,m)';linspace(1/n,1,n)};
    domain = Domain(x);

    verifyEqual(testCase, domain.x, {linspace(1/m,1,m)';linspace(1/n,1,n)})
    verifyEqual(testCase, domain.shape, [m, n])
    verifyEqual(testCase, domain.dimension, 2)
end

function testDomain3dConstructor(testCase)
    l = 256;
    m = 256;
    n = 256;
    x = {linspace(1/l,1,l);linspace(1/m,1,m);linspace(1/n,1,n)};
    domain = Domain(x);

    x3 = zeros(1,1,n);
    x3(1,1,:) = linspace(1/n,1,n);

    verifyEqual(testCase, domain.x, {linspace(1/l,1,l)';linspace(1/m,1,m);x3});
    verifyEqual(testCase, domain.shape, [l, m, n]);
    verifyEqual(testCase, domain.dimension, 3);
end

function testDomain1dReshapeToVector(testCase)
    domain = setup1DDomain();
    f = domain.x{1};

    actual = domain.reshapeToVector(f);
    expected = f;

    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToVector(testCase)
    domain = setup2DDomain();
    f = domain.x{1} + domain.x{2};

    actual = domain.reshapeToVector(f);
    expected = reshape(f,[],1);

    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToVectorLargerHorizontal(testCase)
    domain = setup2DDomain();
    f = domain.x{1} + domain.x{2};
    f = [f, 2*f];

    actual = domain.reshapeToVector(f);
    expected = reshape(f,[],2);

    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToVectorLargerVertical(testCase)
    domain = setup2DDomain();
    f = domain.x{1} + domain.x{2};
    f2 = [f; 2*f];

    actual = domain.reshapeToVector(f2);
    expected = [reshape(f,[],1);2*reshape(f,[],1)];

    verifyEqual(testCase, actual, expected);
end

function testDomain1dReshapeToDomain(testCase)
    domain = setup1DDomain();
    f = domain.x{1};
    fVec = reshape(f,[],1);

    actual = domain.reshapeToDomain(fVec);
    expected = f;

    verifyEqual(testCase, actual, expected);
end

function testDomain1dReshapeToDomainLargerHorizontal(testCase)
    domain = setup1DDomain();
    f = domain.x{1};
    fVec = reshape(f,[],1);
    fVec = [fVec, 2*fVec];

    actual = domain.reshapeToDomain(fVec);
    expected = [f, 2*f];

    verifyEqual(testCase, actual, expected);
end

function testDomain1dReshapeToDomainLargerVertical(testCase)
    domain = setup1DDomain();
    f = domain.x{1};
    fVec = reshape(f,[],1);
    fVec = [fVec; 2*fVec];

    actual = domain.reshapeToDomain(fVec);
    expected = [f; 2*f];

    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToDomain(testCase)
    domain = setup2DDomain();
    f = domain.x{1} + domain.x{2};
    fVec = reshape(f,[],1);

    actual = domain.reshapeToDomain(fVec);
    expected = f;

    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToDomainLargerHorizontal(testCase)
    domain = setup2DDomain();
    f = domain.x{1} + domain.x{2};
%     fVec = domain.reshapeToVector(f);
    fVec = reshape(f,[],1);
    fVec = [fVec, 2*fVec];

    actual = domain.reshapeToDomain(fVec);
    expected = [f, 2*f];

    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToDomainLargerVertical(testCase)
    domain = setup2DDomain();
    f = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
%     fVec = domain.reshapeToVector(f);
    fVec = reshape(f,[],1);
    fVec = [fVec; 2*fVec];

    expected = [f; 2*f];
    actual = domain.reshapeToDomain(fVec);



    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeLargerVertical(testCase)
    domain = setup2DDomain();
    expected = [cos(domain.x{1}) + exp(domain.x{2}); cosh(domain.x{1}+domain.x{2})];

    fVec = domain.reshapeToVector(expected);

    actual = domain.reshapeToDomain(fVec);

    verifyEqual(testCase, actual, expected);
end

function testDomain1dInterpolation(testCase)
    domain = setup1DDomain();
    x = {linspace(0,1,100)};
    f = cos(2*pi*x{1});

    actual = domain.interp(x,f,'spline');
    expected = cos(2*pi*domain.x{1});

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-9, 'RelTol', 1e-6)
end

function testDomain2dInterpolation(testCase)
    domain = setup2DDomain();
    x = {linspace(0,1,100)', linspace(0,1,100)'};
    f = cos(2*pi*x{1}) + cos(2*pi*x{2}');

    actual = domain.interp(x,f,'spline');
    expected = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-6, 'RelTol', 1e-4)
end

function testExtractSurfacesAsCell1D(testCase)
    domain = Domain({linspace(0,1)});
    testSurface = kron([1,2;3,4],domain.x{1});

    expected = {domain.x{1}, 2*domain.x{1}; 3*domain.x{1}, 4*domain.x{1}};
    actual = domain.extractSurfacesAsCell(testSurface);

    verifyEqual(testCase, actual, expected);
end

function testExtractSurfacesAsCell2D(testCase)
    domain = Domain({linspace(0,1),linspace(0,1)});
    y = domain.x{1} + domain.x{2};
    testSurface = kron([1,2;3,4], y);

    expected = {y, 2*y; 3*y, 4*y};
    actual = domain.extractSurfacesAsCell(testSurface);

    verifyEqual(testCase, actual, expected);
end

function testExtractVectorAsCell(testCase)
    domain = Domain({linspace(0,1),linspace(0,1)});
    y = domain.x{1} + domain.x{2};
    Y = domain.reshapeToVector(y);
    testSurface = kron([1,2;3,4], y);
    testVector = domain.reshapeToVector(testSurface);

    expected = {Y, 2*Y; 3*Y, 4*Y};
    actual = domain.extractVectorAsCell(testVector);

    verifyEqual(testCase, actual, expected);
end

%% Functions
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
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});

    dy1 = domain.diff(y, [1, 0]');
    dy2 = domain.diff(y, [0, 1]');

    actual = dy1 / 2 / pi + dy2 / 2 / pi;
    expected = -sin(2*pi*domain.x{1}) - sin(2*pi*domain.x{2});
end

function [actual, expected] = function2dVectorised(domain)
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    y2 = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});

    dy1 = domain.diff([y, y2], [1, 0]');
    dy2 = domain.diff([y, y2], [0, 1]');

    actual = dy1 / 2 / pi + dy2 / 2 / pi;
    expected = cat(3, ...
        -sin(2*pi*domain.x{1})-sin(2*pi*domain.x{2}), ...
        -sin(2*pi*domain.x{1})-sin(2*pi*domain.x{2}));
end

function domain = setup1DDomain()
    n = 256;
    x = {linspace(1/n, 1, n)'};
    domain = Domain(x);
end

function domain = setup2DDomain()
    m = 256;
    n = 256;
    x = {linspace(1/m,1,m)';linspace(1/n,1,n)};
    domain = Domain(x);
end
