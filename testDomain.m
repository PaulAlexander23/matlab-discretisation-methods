function tests = testDomain()
    tests = functiontests(localfunctions);
end

%% Domain
function testDomain1dConstructor(testCase)
    domain = Domain(setup1dX(2^8));
    
    actualX = domain.x;
    actualShape = domain.shape;
    
    expectedX = setup1dX(2^8);
    expectedShape = [2^8, 1];
    
    verifyEqual(testCase, actualX, expectedX)
    verifyEqual(testCase, actualShape, expectedShape)
    verifyEqual(testCase, domain.dimension, 1)
end

function testDomain2dConstructor(testCase)
    domain = Domain(setup2dX(2^8));
    
    actualX = domain.x;
    actualShape = domain.shape;
    
    expectedX = setup2dX(2^8);
    expectedShape = [2^8, 2^8];
    
    verifyEqual(testCase, actualX, expectedX)
    verifyEqual(testCase, actualShape, expectedShape)
    verifyEqual(testCase, domain.dimension, 2)
end

function testDomain1dReshapeToVector(testCase)
    domain = Domain(setup1dX(2^8));
    f = domain.x{1};
    
    actual = domain.reshapeToVector(f);
    expected = f;
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToVector(testCase)
    domain = Domain(setup2dX(2^8));
    f = domain.x{1} + domain.x{2};
    
    actual = domain.reshapeToVector(f);
    expected = reshape(f,[],1);
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToVectorLargerHorizontal(testCase)
    domain = Domain(setup2dX(2^8));
    f = domain.x{1} + domain.x{2};
    f = [f, 2*f];
    
    actual = domain.reshapeToVector(f);
    expected = reshape(f,[],2);
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToVectorLargerVertical(testCase)
    domain = Domain(setup2dX(2^8));
    f = domain.x{1} + domain.x{2};
    f2 = [f; 2*f];
    
    actual = domain.reshapeToVector(f2);
    expected = [reshape(f,[],1);2*reshape(f,[],1)];
    
    verifyEqual(testCase, actual, expected);
end

function testDomain1dReshapeToDomain(testCase)
    domain = Domain(setup1dX(2^8));
    f = domain.x{1};
    fVec = reshape(f,[],1);
    
    actual = domain.reshapeToDomain(fVec);
    expected = f;
    
    verifyEqual(testCase, actual, expected);
end

function testDomain1dReshapeToDomainLargerHorizontal(testCase)
    domain = Domain(setup1dX(2^8));
    f = domain.x{1};
    fVec = reshape(f,[],1);
    fVec = [fVec, 2*fVec];
    
    actual = domain.reshapeToDomain(fVec);
    expected = [f, 2*f];
    
    verifyEqual(testCase, actual, expected);
end

function testDomain1dReshapeToDomainLargerVertical(testCase)
    domain = Domain(setup1dX(2^8));
    f = domain.x{1};
    fVec = reshape(f,[],1);
    fVec = [fVec; 2*fVec];
    
    actual = domain.reshapeToDomain(fVec);
    expected = [f; 2*f];
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToDomain(testCase)
    domain = Domain(setup2dX(2^8));
    f = domain.x{1} + domain.x{2};
    fVec = reshape(f,[],1);
    
    actual = domain.reshapeToDomain(fVec);
    expected = f;
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToDomainLargerHorizontal(testCase)
    domain = Domain(setup2dX(2^8));
    f = domain.x{1} + domain.x{2};
%     fVec = domain.reshapeToVector(f);
    fVec = reshape(f,[],1);
    fVec = [fVec, 2*fVec];
    
    actual = domain.reshapeToDomain(fVec);
    expected = [f, 2*f];
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeToDomainLargerVertical(testCase)
    domain = Domain(setup2dX(2^8));
    f = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
%     fVec = domain.reshapeToVector(f);
    fVec = reshape(f,[],1);
    fVec = [fVec; 2*fVec];
    
    expected = [f; 2*f];
    actual = domain.reshapeToDomain(fVec);
    
    
    
    verifyEqual(testCase, actual, expected);
end

function testDomain2dReshapeLargerVertical(testCase)
    domain = Domain(setup2dX(2^8));
    expected = [cos(domain.x{1}) + exp(domain.x{2}); cosh(domain.x{1}+domain.x{2})];
    
    fVec = domain.reshapeToVector(expected);
    
    actual = domain.reshapeToDomain(fVec);
    
    verifyEqual(testCase, actual, expected);
end

function testDomain1dInterpolation(testCase)
    domain = Domain(setup1dX(2^8));
    x = {linspace(0,1,100)};
    f = cos(2*pi*x{1});
    
    actual = domain.interp(x,f);
    expected = cos(2*pi*domain.x{1});
    
    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-9, 'RelTol', 1e-6)
end

function testDomain2dInterpolation(testCase)
    domain = Domain(setup2dX(2^8));
    x = {linspace(0,1,100)', linspace(0,1,100)'};
    f = cos(2*pi*x{1}) + cos(2*pi*x{2}');
    
    actual = domain.interp(x,f);
    expected = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-6, 'RelTol', 1e-4)
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

function x = setup1dX(n)
    x = {linspace(1/n, 1, n)'};
end

function x = setup2dX(m,n)
    if nargin < 2
        n = m;
    end
    x = {linspace(1/m,1,m)';linspace(1/n,1,n)};
end