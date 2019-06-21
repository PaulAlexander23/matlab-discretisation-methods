function tests = tests()
    tests = functiontests(localfunctions);
end

function testEvaluatingFunction1dFiniteDifference(testCase)
    n = 2^8;
    x = linspace(1/n,1,n)';
    
    problemDeg = [1, 2];
    
    D = init_fd({x}, problemDeg, 2);
    method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
    [actual, expected] = problem1d(x,method);
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction1dPseudoSpectral(testCase)
    n = 2^8;
    x = linspace(1/n,1,n)';
    method = @diff_ps;
    [actual, expected] = problem1d(x,method);
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
end

function testEvaluatingFunction2dFiniteDifference(testCase)
    n = 2^8;
    x = linspace(1/n,1,n)';
    
    problemDeg = [1,0;0,1]';
    
    D = init_fd({x}, problemDeg, 2);
    method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
    [actual, expected] = problem2d(x,method);
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction2dPseudoSpectral(testCase)
    n = 2^8;
    x = linspace(1/n,1,n)';
    method = @diff_ps;
    [actual, expected] = problem2d(x,method);
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
end

function testEvaluatingFunction2dDivFiniteDifference(testCase)
    n = 2^8;
    x = linspace(1/n,1,n)';
    
    problemDeg = cat(3,[1,0;0,0],[0,0;0,1]);
    
    D = init_fd({x}, problemDeg, 2);
    method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
    [actual, expected] = problem2dDiv(x,method);
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction2dDivPseudoSpectral(testCase)
    n = 2^8;
    x = linspace(1/n,1,n)';
    method = @diff_ps;
    [actual, expected] = problem2dDiv(x,method);
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
end

function problemMethod(testCase, problem, problemDeg, expectedFunction, diffMethod)
    n = 2^8;
    x = linspace(1/n,1,n)';
    
    y = cos(2*pi*x);
    
    if diffMethod == "diff_fd"
        D = init_fd({x}, problemDeg, 2);
        method = @(x,y,deg) diff_fd(x,y,deg,D,problemDeg);
    elseif diffMethod == "diff_ps"
        method = @diff_ps;
    end
    
    actual = feval(problem,{x},y,method);
    expected = expectedFunction(x);
    
    if diffMethod == "diff_fd"
        verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
    elseif diffMethod == "diff_ps"
        verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
    end
end

function [actual, expected] = function1d(x, method)
    deg = [1, 2];
    y = cos(2 * pi * x{1});
    dy = method(x, y, deg);
    actual = dy{2}/(4*pi^2) + dy{1}/(2*pi) + y;
    expected = - sin(2 * pi * x);
end

function [actual, expected] = function2d(x,method)
    deg = [1,0;0,1]';
    y = cos(2 * pi * x{1}) + cos(2 * pi * x{2});
    dy = method(x,y,deg);
    actual = dy{1}/2/pi + dy{2}/pi;
    expected = -sin(2 * pi * x{1}) - 2 * sin(2 * pi * x{2});
end

function [actual, expected] = function2dDiv(x,method)
    divdeg = cat(3,[1,0;0,0],[0,0;0,1]);
    y = cos(2 * pi * x{1}) + cos(2 * pi * x{2});
    dy = method(x,y,divdeg);
    actual = dy{1} + dy{2};
    expected = cat(3, - 2 * pi * sin(2 * pi * x{1}) * ones(size(x{2})), ...
        - 2 * pi * sin(2 * pi * x{2}) * ones(size(x{1})));
end