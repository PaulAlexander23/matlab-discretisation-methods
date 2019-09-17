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
    Y = domain.fft(cos(2*pi*domain.x{1}));

    degree = 1;
    
    actual = domain.ifft(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*domain.x{1});
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-14)
end

function test2dPseudoSpectralDiff(testCase)
    domain = PSDomain(setup2dX(2^8, 2^7));
    Y = domain.fft(cos(2*pi*domain.x{1}) + 0 * domain.x{2}');

    actual = domain.ifft(domain.diff(Y, [1, 0]'));
    expected = -2*pi*sin(2*pi*domain.x{1}) + 0 * domain.x{2}';
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-14)
end

function test1dFiniteDifferenceGetDiffMatrix(testCase)
    domain = FDDomain(setup1dX(2^8), 1, 2);
    degree = 1;
    
    actual = domain.getDiffMatrix(degree);

    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize);
end

function testPseudoSpectralScaling1d(testCase)
    n = 2^8;
   domain = PSDomain(setup1dX(n));
   amplitude = 2;
   y = amplitude * cos(2 * pi * domain.x{1});
   fy = domain.fft(y);
   verifyEqual(testCase, max(real(fy)), amplitude, 'RelTol', 1e-3);
end

function testPseudoSpectralScaling2d(testCase)
   n = 2^8;
   domain = PSDomain(setup2dX(n));
   amplitude = 2;
   y = amplitude * cos(2 * pi * (domain.x{1} + domain.x{2}'));
   fy = domain.fft(y);
   verifyEqual(testCase, max(max(real(fy))), amplitude, 'RelTol', 1.3e-3);
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
    domain = PSDomain({linspace(1/2^7,1,2^7)';linspace(1/2^8,1,2^8)'});
    Y = domain.fft(cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}'));
    degree = [1, 0]';
    
    actual = domain.ifft(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*domain.x{1}) .* ones(size(domain.x{2}'));
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14)
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
    
    verifySize(testCase, actual, [2^8, 2^8, 2])
    verifyEqual(testCase,actual,expected,'RelTol',1e-3,'AbsTol',1e-4)
end

function testZeropad1D(testCase)
    N = 32;
    domain = PSDomain(setup1dX(N));
    
    f = domain.x{1};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    
    expected = zeros(48,1);
    expected(1:N/2) = f(1:N/2);
    expected(ratio*N-N/2+1:ratio*N) = f(1+N/2:end);
    
    verifyEqual(testCase, actual, expected);
end

function testZeropad2D(testCase)
    N = 32;
    domain = PSDomain(setup2dX(N));
    
    f = domain.x{1} + domain.x{2}';
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    
    expected = zeros(48);
    expected(1:N/2,1:N/2) = f(1:N/2,1:N/2);
    expected(1:N/2,ratio*N-N/2+1:ratio*N) = f(1:N/2,1+N/2:end);
    expected(ratio*N-N/2+1:ratio*N,1:N/2) = f(1+N/2:end,1:N/2);
    expected(ratio*N-N/2+1:ratio*N,ratio*N-N/2+1:ratio*N) = f(1+N/2:end,1+N/2:end);
    
    verifyEqual(testCase, actual, expected);
end

function testTruncation1D(testCase)
    N = 3*32;
    domain = PSDomain(setup1dX(N));
    
    f = domain.x{1};
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    
    expected = f([1:N*1/3, N*2/3+1:N]);
    
    verifyEqual(testCase, actual, expected);
end

function testTruncation2D(testCase)
    N = 3*32;
    domain = PSDomain(setup2dX(N));
    
    f = domain.x{1} + domain.x{2}';
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    
    expected = f([1:N*1/3, N*2/3+1:N], [1:N*1/3, N*2/3+1:N]);
    
    verifyEqual(testCase, actual, expected);
end

function testConvolution1d(testCase)
    N = 2^7;
    domain = PSDomain(setup1dX(N));
    y = cos(2*pi*domain.x{1});
    fy = domain.fft(y);
    actual = domain.multiply(fy, fy);
    expected = zeros(N,1);
    expected(1) = 1;
    expected(3) = 0.5;
    expected(N-1) = 0.5;
    
    verifyEqual(testCase, real(actual), expected, 'AbsTol', 1e-16, 'RelTol', 5e-3);
end

function testConvDealiasing(testCase)
    N = 2^5; M = 2^6;
    coarseDomain = PSDomain(setup1dX(N));
    fineDomain = PSDomain(setup1dX(M));
    u = sin(N/2 * pi * fineDomain.x{1});
    v = cos((N/2 + 2) * pi * fineDomain.x{1});
    uhat = fineDomain.fft(u);
    vhat = fineDomain.fft(v);

    actual = coarseDomain.multiply(uhat([1:N/2,M-N/2+1:M]), vhat([1:N/2,M-N/2+1:M]));
    
    expected = fineDomain.fft(u .* v);

%     plot(real(actual))
%     hold on
%     plot(real(expected([1:N/2,M-N/2+1:M])))
    
    verifyEqual(testCase, actual, expected([1:N/2,M-N/2+1:M]), 'AbsTol', 1e-13);
end

function testConvDealiasingHigherPowers(testCase)
    N = 2^7; M = 2^8;
    coarseDomain = PSDomain(setup1dX(N));
    fineDomain = PSDomain(setup1dX(M));
    u = sin(ceil(N/3) * pi * fineDomain.x{1});
    v = cos((ceil(N/3) + 2) * pi * fineDomain.x{1});
    uhat = fineDomain.fft(u);
    vhat = fineDomain.fft(v);
    
    actual = coarseDomain.multiply(uhat([1:N/2,M-N/2+1:M]), vhat([1:N/2,M-N/2+1:M]), [2, 1]);
    expected = fineDomain.fft(u.^2 .* v);

%     figure;
%     hold on;
%     plot(log10(abs(actual - expected([1:N/2,M-N/2+1:M]))));
%     plot(log10(abs(actual - expected([1:N/2,M-N/2+1:M]))./abs(actual)));
    
    verifyEqual(testCase, actual, expected([1:N/2,M-N/2+1:M]), 'AbsTol', 1e-2, 'RelTol', 1e-1);
end

function testConvDealiasingNegativePowers(testCase)
    N = 2^7; M = 2^8;
    coarseDomain = PSDomain(setup1dX(N));
    fineDomain = PSDomain(setup1dX(M));
    u = sin(ceil(N/3) * pi * fineDomain.x{1}) + 2;
    v = cos((ceil(N/3) + 2) * pi * fineDomain.x{1});
    uhat = fineDomain.fft(u);
    vhat = fineDomain.fft(v);

    actual = coarseDomain.multiply(uhat([1:N/2,M-N/2+1:M]), vhat([1:N/2,M-N/2+1:M]), [-1, 1]);
    expected = fineDomain.fft((fineDomain.ifft(uhat)).^-1 .* (fineDomain.ifft(vhat)));
    
%     figure;
%     plot(real(actual))
%     hold on
%     plot(real(expected([1:N/2,M-N/2+1:M])))
%     figure;
%     hold on;
%     plot(log10(abs(actual - expected([1:N/2,M-N/2+1:M]))));
%     plot(log10(abs(actual - expected([1:N/2,M-N/2+1:M]))./abs(actual)));

    verifyEqual(testCase, actual, expected([1:N/2,M-N/2+1:M]), 'AbsTol', 1e-2, 'RelTol', 1e-1);
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

function x = setup2dX(m,n)
    if nargin < 2
        n = m;
    end
    x = {linspace(1/m,1,m)';linspace(1/n,1,n)'};
end
