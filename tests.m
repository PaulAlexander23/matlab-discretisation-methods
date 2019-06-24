function tests = tests()
    tests = functiontests(localfunctions);
end

function testDomain1d(testCase)
    actual = setup1dX(2^8);
    
    domain = Domain(actual);
    expected = domain.x;
    
    verifyEqual(testCase, actual, expected)
end

function testDomain2d(testCase)
    actual = setup2dX(2^8);
    
    domain = Domain(actual);
    expected = domain.x;
    
    verifyEqual(testCase, actual, expected)
end

function test1dFiniteDifference(testCase)
    x = setup1dX(2^8);
    
    domain = FDDomain(x, 1, 2);
    
    Y = cos(2*pi*domain.x{1});
    
    degree = 1;
    
    actual = cell2mat(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*x{1});
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test1dPseudoSpectral(testCase)
    x = setup1dX(2^8);
    
    domain = PSDomain(x);
    
    Y = cos(2*pi*domain.x{1});
    
    degree = 1;
    
    actual = cell2mat(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*x{1});
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-14)
end

function testEvaluatingFunction1dFiniteDifference(testCase)
    x = setup1dX(2^8);
    
    domain = FDDomain(x, [1, 2], 2);
    
    method = @(~, y, deg) domain.diff(y, deg);
    
    [actual, expected] = function1d(domain.x,method);
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction1dPseudoSpectral(testCase)
    x = setup1dX(2^8);
    
    domain = PSDomain(x);
    
    method = @(~, y, deg) domain.diff(y, deg);
    
    [actual, expected] = function1d(domain.x,method);
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
end

function test2dFiniteDifference(testCase)
    x = setup2dX(2^8);
    
    domain = FDDomain(x, [1, 0]', 2);
    
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');
    Y = reshape(Y,[1, numel(Y)]);
    
    degree = [1, 0]';
    
    actual = cell2mat(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*x{1}) .* ones(size(x{2}'));
    
    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test2dPseudoSpectral(testCase)
    x = setup2dX(2^8);
    
    domain = PDDomain(x, [1, 0]', 2);
    
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');
    Y = reshape(Y,[1, numel(Y)]);
    
    degree = [1, 0]';
    
    actual = cell2mat(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*x{1}) .* ones(size(x{2}'));
    
    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function testEvaluatingFunction2dFiniteDifference(testCase)
    x = setup2dX(2^8);
    
    problemDeg = [1,0;0,1]';
    
    domain = FDDomain(x, problemDeg, 2);
    
    method = @(~, y, deg) domain.diff(y, deg);
    
    [actual, expected] = function2d(domain.x,method);
    
    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction2dPseudoSpectral(testCase)
    x = setup2dX(2^8);
    
    domain = PSDomain(x);
    
    method = @(~, y, deg) domain.diff(y, deg);
    
    [actual, expected] = function2d(domain.x,method);
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
end

function [actual, expected] = function1d(x, method)
    y = cos(2 * pi * x{1});
    dy = method(x, y, [1, 2]);
    
    actual = dy{2}/(4*pi^2) + dy{1}/(2*pi) + y;
    expected = - sin(2 * pi * x{1});
end

function [actual, expected] = function2d(x,method)
    y = cos(2 * pi * x{1}) + cos(2 * pi * x{2}');
    y = reshape(y,[1, numel(y)]);
    
    dy1 = cell2mat(method(x, y, [1, 0]'));
    dy2 = cell2mat(method(x, y, [0, 1]'));
    
    actual = dy1/2/pi + dy2/2/pi;
    expected = -sin(2 * pi * x{1}) - sin(2 * pi * x{2}');
end

function x = setup1dX(n)
    x = {linspace(1/n,1,n)'};
end

function x = setup2dX(n)
    x = {linspace(1/n,1,n)';linspace(1/n,1,n)'};
end