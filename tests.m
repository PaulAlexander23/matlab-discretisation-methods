function tests = tests()
    tests = functiontests(localfunctions);
end

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

function test1dPseudoSpectralConstructor(testCase)
    domain = PSDomain(setup1dX(2^8));
    
    verifyEqual(testCase, domain.length, 1);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        [0:2^7 - 1, 0, 1 - 2^7:-1]' * 2*pi);
end

function test1dFiniteDifferenceDiff(testCase)
    domain = FDDomain(setup1dX(2^8), 1, 2);
    Y = cos(2*pi*domain.x{1});
    degree = 1;
    
    actual = domain.diff(Y, degree);

    expected = -2*pi*sin(2*pi*domain.x{1});
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function test1dPseudoSpectralDiff(testCase)
    domain = PSDomain(setup1dX(2^8));
    Y = fft(cos(2*pi*domain.x{1}));
    degree = 1;
    
    actual = ifft(domain.diff(Y, degree));

    expected = -2*pi*sin(2*pi*domain.x{1});
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-14)
end

function test1dFiniteDifferenceGetDiffMatrix(testCase)
    domain = FDDomain(setup1dX(2^8), 1, 2);
    degree = 1;
    
    actual = domain.getDiffMatrix(degree);

    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize);
end

function testEvaluatingFunction1dFiniteDifference(testCase)
    domain = FDDomain(setup1dX(2^8), [1, 2], 2);
    
    [actual, expected] = function1d(domain);
    
    verifySize(testCase, actual, [2^8, 1])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testEvaluatingFunction1dVectorisedFiniteDifference(testCase)
    domain = FDDomain(setup1dX(2^8), [1, 2], 2);
    
    [actual, expected] = function1dVectorised(domain);
    
    verifySize(testCase, actual, [2^8, 2])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

% function testEvaluatingFunction1dPseudoSpectral(testCase)
%     x = setup1dX(2^8);
%     
%     domain = PSDomain(x);
%     
%     [actual, expected] = function1d(domain);
%     
%     verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-15)
% end

function test2dFiniteDifference(testCase)
    domain = FDDomain(setup2dX(2^8), [1, 0]', 2);
    Y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');
    degree = [1, 0]';
    
    actual = domain.diff(Y, degree);

    expected = -2*pi*sin(2*pi*domain.x{1}) .* ones(size(domain.x{2}'));
    
    verifySize(testCase, actual, [2^8, 2^8])
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-3, 'AbsTol', 1e-4)
end

function test2dPseudoSpectral(testCase)
    domain = PSDomain(setup2dX(2^8));
    Y = fftn(cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    degree = [1, 0]';
    
    actual = ifftn(domain.diff(Y, degree));

    expected = -2*pi*sin(2*pi*domain.x{1}) .* ones(size(domain.x{2}'));
    
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

% function testEvaluatingFunction2dPseudoSpectral(testCase)
%     x = setup2dX(2^8);
%     
%     domain = PSDomain(x);
%     
%     [actual, expected] = function2d(domain);
%     
%     verifyEqual(testCase,actual,expected,'RelTol',1e-12,'AbsTol',1e-13)
% end

function testEvaluatingFunction2dVectorisedFiniteDifference(testCase)
    problemDeg = [1,0;0,1]';
    domain = FDDomain(setup2dX(2^8), problemDeg, 2);
    
    [actual, expected] = function2dVectorised(domain);
    
    verifySize(testCase, actual, [2^8, 2^8, 2])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testFiniteDifferenceFbenney2d(testCase)
    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2]';
    domain = FDDomain(setup2dX(2^8), diffDegrees, 4);
    y = cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}');
    params = [1, 7/8*pi, 1, 0.01];
    
    actual = fbenney2d(domain, y, params);

    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize)
end

function testPseudoSpectralFbenney2d(testCase)
    domain = PSDomain(setup2dX(2^8));
    y = domain.fftn(cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    params = [1, 7/8*pi, 1, 0.01];
    
    actual = fbenney2d(domain, y, params);

    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize)
end

function testConvDealiasing(testCase)
    N = 2^5; M = 2^6;
    coarseDomain = PSDomain(setup1dX(N));
    fineDomain = PSDomain(setup1dX(M));
    u = sin(N/2 * pi * fineDomain.x{1});
    v = cos((N/2 + 2) * pi * fineDomain.x{1});
    uhat = fineDomain.fftn(u);
    vhat = fineDomain.fftn(v);

    actual = coarseDomain.multiply(uhat([1:N/2,M-N/2+1:M]), vhat([1:N/2,M-N/2+1:M]));
    
    expected = fineDomain.fftn(fineDomain.ifftn(uhat) .* fineDomain.ifftn(vhat));

    plot(real(actual))
    hold on
    plot(real(expected([1:N/2,M-N/2+1:M])))
    
    verifyEqual(testCase, actual, expected([1:N/2,M-N/2+1:M]), 'AbsTol', 1e-13);
end

function testConvDealiasingHigherPowers(testCase)
    N = 2^5; M = 2^6;
    coarseDomain = PSDomain(setup1dX(N));
    fineDomain = PSDomain(setup1dX(M));
    u = sin(ceil(N/3) * pi * fineDomain.x{1});
    v = cos((ceil(N/3) + 2) * pi * fineDomain.x{1});
    uhat = fineDomain.fftn(u);
    vhat = fineDomain.fftn(v);

    actual = coarseDomain.multiply(uhat([1:N/2,M-N/2+1:M]), vhat([1:N/2,M-N/2+1:M]), [2, 1]);
    expected = fineDomain.fftn(fineDomain.ifftn(uhat).^2 .* fineDomain.ifftn(vhat));
    
    plot(real(actual))
    hold on
    plot(real(expected([1:N/2,M-N/2+1:M])))

    verifyEqual(testCase, actual, expected([1:N/2,M-N/2+1:M]), 'AbsTol', 1e-2, 'RelTol', 1e-2);
end

function testConvDealiasingNegativePowers(testCase)
    N = 2^5; M = 2^6;
    coarseDomain = PSDomain(setup1dX(N));
    fineDomain = PSDomain(setup1dX(M));
    u = sin(ceil(N/3) * pi * fineDomain.x{1}) + 2;
    v = cos((ceil(N/3) + 2) * pi * fineDomain.x{1});
    uhat = fft(u)/M;
    vhat = fft(v)/M;

    actual = coarseDomain.multiply(uhat([1:N/2,M-N/2+1:M]), vhat([1:N/2,M-N/2+1:M]), [-1, 1]);
    expected = fineDomain.fftn((fineDomain.ifftn(uhat)).^-1 .* (fineDomain.ifftn(vhat)));
    
    plot(real(actual))
    hold on
    plot(real(expected([1:N/2,M-N/2+1:M])))

    verifyEqual(testCase, actual, expected([1:N/2,M-N/2+1:M]), 'AbsTol', 2e-2, 'RelTol', 1e-2);
end

function testFiniteDifferenceDealiasingOnWIBL1(testCase)
    domain = FDDomain(setup2dX(2^6),[1,0;0,1;2,0;0,2]',4);
    y = 1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
    f = 2/3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}');
    Y = [y;f];
    params = [1, pi/4, 1, 0.01];
    
    actual = fwibl1(domain, Y, params);

    verifyTrue(testCase, all(max(actual)<1e5))
end

function testPseudoSpectralDealiasingOnWIBL1(testCase)
    domain = PSDomain(setup2dX(2^6));
    y = domain.fftn(1 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}'));
    f = domain.fftn(2/3 + 0.1 * cos(4*pi*domain.x{1}) + 0.1 * cos(4*pi*domain.x{2}'));
    Y = [y;f];
    params = [1, pi/4, 1, 0.01];
    
    actual = fwibl1(domain, Y, params);
    
    verifyTrue(testCase, all(max(actual)<1e5))
end

function [actual, expected] = function1d(domain)
    y = cos(2 * pi * domain.x{1});
    dy = domain.diff(y, 1);
    d2y = domain.diff(y, 2);
    actual = d2y/(4*pi^2) + dy/(2*pi) + y;
    expected = - sin(2 * pi * domain.x{1});
end

function [actual, expected] = function1dVectorised(domain)
    y = cos(2 * pi * domain.x{1});
    y2 = cos(2 * pi * domain.x{1});
    dy = domain.diff([y, y2], 1);
    d2y = domain.diff([y, y2], 2);
    actual = d2y/(4*pi^2) + dy/(2*pi) + y;
    expected = cat(2, ...
        - sin(2 * pi * domain.x{1}), ...
        - sin(2 * pi * domain.x{1}));
end

function [actual, expected] = function2d(domain)
    y = cos(2 * pi * domain.x{1}) + cos(2 * pi * domain.x{2}');
    
    dy1 = domain.diff(y, [1, 0]');
    dy2 = domain.diff(y, [0, 1]');
    
    actual = dy1/2/pi + dy2/2/pi;
    expected = -sin(2 * pi * domain.x{1}) - sin(2 * pi * domain.x{2}');
end

function [actual, expected] = function2dVectorised(domain)
    y = cos(2 * pi * domain.x{1}) + cos(2 * pi * domain.x{2}');
    y2 = cos(2 * pi * domain.x{1}) + cos(2 * pi * domain.x{2}');
    
    dy1 = domain.diff([y, y2], [1, 0]');
    dy2 = domain.diff([y, y2], [0, 1]');
    
    actual = dy1/2/pi + dy2/2/pi;
    expected = cat(3, ...
        -sin(2 * pi * domain.x{1}) - sin(2 * pi * domain.x{2}'), ...
        -sin(2 * pi * domain.x{1}) - sin(2 * pi * domain.x{2}'));
end

function x = setup1dX(n)
    x = {linspace(1/n,1,n)'};
end

function x = setup2dX(n)
    x = {linspace(1/n,1,n)';linspace(1/n,1,n)'};
end
