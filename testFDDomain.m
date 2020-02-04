function tests = testFDDomain()
    tests = functiontests(localfunctions);
end

%% Finite Difference
function test1dFiniteDifferenceDiff(testCase)
    domain = FDDomain(setup1dX(2^8), 1, 2);
    Y = cos(2*pi*domain.x{1});
    degree = 1;
    
    actual = domain.diff(Y, degree);
    
    expected = -2*pi*sin(2*pi*domain.x{1});
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test1dFiniteDifferenceDiffPair(testCase)
    domain = FDDomain(setup1dX(2^8), 1, 2);
    Y = cos(2*pi*domain.x{1});
    Y = [Y;2*Y];
    degree = 1;
    
    actual = domain.diff(Y, degree);
    
    expected = [-2*pi*sin(2*pi*domain.x{1});-4*pi*sin(2*pi*domain.x{1})];
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test1dFiniteDifferenceGetDiffMatrix(testCase)
    domain = FDDomain(setup1dX(2^8), 1, 2);
    degree = 1;
    
    actual = domain.diffMat(degree);
    
    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize);
end

function test2dFiniteDifferenceDiff(testCase)
    domain = FDDomain(setup2dX(2^8), [1,0;0,1]', 2);
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    
    expected = -2*pi*sin(2*pi*domain.x{1}) + 0*domain.x{2};
    actual = domain.diff(Y, [1,0]');
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
    
    expected = -2*pi*sin(2*pi*domain.x{2}) + 0*domain.x{1};
    actual = domain.diff(Y, [0,1]');
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test2dFiniteDifferenceDiffPair(testCase)
    domain = FDDomain(setup2dX(2^8), [1,0;0,1]', 2);
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    Y = [Y; 2*Y];
    
    expected = [-2*pi*sin(2*pi*domain.x{1}) + 0*domain.x{2};
        -4*pi*sin(2*pi*domain.x{1}) + 0*domain.x{2}];
    actual = domain.diff(Y, [1,0]');
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)

    expected = [-2*pi*sin(2*pi*domain.x{2}) + 0*domain.x{1}; ...
        -4*pi*sin(2*pi*domain.x{2}) + 0*domain.x{1}];
    actual = domain.diff(Y, [0,1]');
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test2dFiniteDifferenceDiffVector(testCase)
    domain = FDDomain(setup2dX(2^8), [1,0;0,1]', 2);
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    Y = [Y, 2*Y];
    
    expected = [-2*pi*sin(2*pi*domain.x{1}) + 0*domain.x{2},...
        -4*pi*sin(2*pi*domain.x{1}) + 0*domain.x{2}];
    actual = domain.diff(Y, [1,0]');
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
    
    expected = [-2*pi*sin(2*pi*domain.x{2}) + 0*domain.x{1},...
        -4*pi*sin(2*pi*domain.x{2}) + 0*domain.x{1}];
    actual = domain.diff(Y, [0,1]');
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test2dFiniteDifferenceGetDiffMatrix(testCase)
    degree =  [1, 0]';
    domain = FDDomain(setup2dX(2^8), degree, 2);
    
    actual = domain.diffMat(degree);
    
    expectedSize = [2^16, 2^16];
    
    verifySize(testCase, actual, expectedSize);
end

function testEvaluatingFunction1dFiniteDifference(testCase)
    domain = FDDomain(setup1dX(2^8), [1, 2], 2);
    
    [actual, expected] = function1d(domain);
    
    verifySize(testCase, actual, [2^8, 1])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction1dFiniteDifferenceForwardDifference(testCase)
    domain = FDDomain(setup1dX(2^8), [1, 2], 2, 'forward');
    
    [actual, expected] = function1d(domain);
    
    verifySize(testCase, actual, [2^8, 1])
    verifyEqual(testCase,actual,expected,'RelTol',1e-2,'AbsTol',1e-3)
end

function testEvaluatingFunction1dVectorisedFiniteDifference(testCase)
    domain = FDDomain(setup1dX(2^8), [1, 2], 2);
    
    [actual, expected] = function1dVectorised(domain);
    
    verifySize(testCase, actual, [2^8, 2])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test2dFiniteDifference(testCase)
    domain = FDDomain(setup2dX(2^8), [1, 0]', 2);
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    degree = [1, 0]';
    
    actual = domain.diff(Y, degree);
    expected = -2*pi*sin(2*pi*domain.x{1}) .* ones(size(domain.x{2}));
    
    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function testEvaluatingFunction2dFiniteDifference(testCase)
    problemDeg = [1,0;0,1]';
    domain = FDDomain(setup2dX(2^8), problemDeg, 2);
    
    [actual, expected] = function2d(domain);
    
    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction2dVectorisedFiniteDifference(testCase)
    problemDeg = [1,0;0,1]';
    domain = FDDomain(setup2dX(2^8), problemDeg, 2);
    
    [actual, expected] = function2dVectorised(domain);
    
    verifySize(testCase, actual, [2^8, 2^8 * 2])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction2dPairFiniteDifference(testCase)
    problemDeg = [1,0;0,1]';
    domain = FDDomain(setup2dX(2^8), problemDeg, 2);
    
    [actual, expected] = function2dPair(domain);
    
    verifySize(testCase, actual, [2^8 * 2, 2^8])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
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
    expected = cat(2, ...
        -sin(2*pi*domain.x{1})-sin(2*pi*domain.x{2}), ...
        -sin(2*pi*domain.x{1})-sin(2*pi*domain.x{2}));
end

function [actual, expected] = function2dPair(domain)
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    z = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2});
    u = [y;z];
    
    dux = domain.diff(u, [1, 0]');
    duy = domain.diff(u, [0, 1]');

    actual = dux / 2 / pi + duy / 2 / pi;
    expected = cat(1, ...
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