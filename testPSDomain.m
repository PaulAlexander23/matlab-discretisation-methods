function tests = testPSDomain()
    tests = functiontests(localfunctions);
end

%% Pseudo Spectral
function test1dPseudoSpectralConstructorComplex(testCase)
    domain = PSDomain(setup1dX(2^8), true, true, 1e-5);
    
    verifyEqual(testCase, domain.length, 1);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        [0:2^7 - 1, 0, 1 - 2^7:-1]' * 2*pi);
    verifyEqual(testCase, domain.antialiasing, true);
    verifyEqual(testCase, domain.complex, true);
    verifyEqual(testCase, domain.suppression, 1e-5);
end

function test1dPseudoSpectralConstructorReal(testCase)
    domain = PSDomain(setup1dX(2^8), true, false, 1e-5);
    
    verifyEqual(testCase, domain.length, 1);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        (0:2^7 - 1)' * 2*pi);
    verifyEqual(testCase, domain.antialiasing, true);
    verifyEqual(testCase, domain.complex, false);
    verifyEqual(testCase, domain.suppression, 1e-5);
end

function test1dPseudoSpectralfft(testCase)
    N = 2^4;
    domain = PSDomain(setup1dX(N));
    y = cos(2*pi*(domain.x{1} - domain.x{1}(1)));
    
    actual = domain.fft(y);
    expected = [0; 1; zeros(N-3, 1); 1];
    
    verifyEqual(testCase, actual, expected, 'RelTol', eps, 'AbsTol', eps);
end

function test1dPseudoSpectralfftSmallestWavenumber(testCase)
    N = 2^4;
    domain = PSDomain(setup1dX(N));
    y = cos((N/2-1)*2*pi*(domain.x{1} - domain.x{1}(1)));
    
    actual = domain.fft(y);
    expected = [0; zeros(N/2-2, 1); 1; 0; 1; zeros(N/2-2, 1)];
    
    verifyEqual(testCase, actual, expected, 'RelTol', eps, 'AbsTol', 1e-15);
end

function test1dPseudoSpectralfftReal(testCase)
    N = 2^4;
    domain = PSDomain(setup1dX(N), false, false);
    y = cos(2*pi*(domain.x{1} - domain.x{1}(1)));
    
    actual = domain.fft(y);
    expected = [0; 1; zeros(N/2-2, 1)];
    
    verifyEqual(testCase, actual, expected, 'RelTol', eps, 'AbsTol', eps);
end

function test1dPseudoSpectralfftRealSmallestWavenumber(testCase)
    N = 2^4;
    domain = PSDomain(setup1dX(N), false, false);
    y = cos((N/2-1)*2*pi*(domain.x{1} - domain.x{1}(1)));
    
    actual = domain.fft(y);
    expected = [0; zeros(N/2-2, 1); 1];
    
    verifyEqual(testCase, actual, expected, 'RelTol', eps, 'AbsTol', 1e-15);
end

function test1dPseudoSpectralifft(testCase)
    N = 2^4;
    domain = PSDomain(setup1dX(N));
    
    actual = domain.ifft([0; 1; zeros(N-3, 1); 1]);
    expected = cos(2*pi*(domain.x{1} - domain.x{1}(1)));
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', eps);
end

function test1dPseudoSpectralifftReal(testCase)
    N = 2^4;
    domain = PSDomain(setup1dX(N), false, false);
    
    actual = domain.ifft([0; 1; zeros(N/2-2, 1)]);
    expected = cos(2*pi*(domain.x{1} - domain.x{1}(1)));
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', eps);
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
    Y = domain.fft(cos(2*pi*domain.x{1}) + 0 * domain.x{2});
    
    actual = domain.ifft(domain.diff(Y, [1, 0]'));
    expected = -2*pi*sin(2*pi*domain.x{1}) + 0 * domain.x{2};
    
    verifyEqual(testCase,actual,expected,'RelTol',1e-14,'AbsTol',1e-14)
end

function test1dPseudoSpectralGetDiffMatrix(testCase)
    domain = PSDomain(setup1dX(2^8));
    degree = 1;
    
    actual = domain.diffMat(degree);
    
    expectedSize = [2^8, 2^8];
    
    verifySize(testCase, actual, expectedSize);
end

function test2dPseudoSpectralGetDiffMatrix(testCase)
    domain = PSDomain(setup2dX(2^8));
    degree = [1, 0]';
    
    actual = domain.diffMat(degree);
    
    expectedSize = [2^16, 2^16];

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
    y = amplitude * cos(2 * pi * (domain.x{1} + domain.x{2}));
    fy = domain.fft(y);
    verifyEqual(testCase, max(max(real(fy))), amplitude, 'RelTol', 1.3e-3);
end

function test2dPseudoSpectral(testCase)
    domain = PSDomain({linspace(1/2^7,1,2^7)';linspace(1/2^8,1,2^8)});
    Y = domain.fft(cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    degree = [1, 0]';
    
    actual = domain.ifft(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*domain.x{1}) .* ones(size(domain.x{2}));
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14)
end

function test2dPseudoSpectralReal(testCase)
    domain = PSDomain({linspace(1/2^7,1,2^7)';linspace(1/2^8,1,2^8)'}, false, false);
    Y = domain.fft(cos(2*pi*domain.x{1}) + cos(2*pi*domain.x{2}));
    degree = [1, 0]';
    
    actual = domain.ifft(domain.diff(Y, degree));
    expected = -2*pi*sin(2*pi*domain.x{1}) .* ones(size(domain.x{2}));
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14)
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
    
    f = domain.x{1} + domain.x{2};
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
    
    f = domain.x{1} + domain.x{2};
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    
    expected = f([1:N*1/3, N*2/3+1:N], [1:N*1/3, N*2/3+1:N]);
    
    verifyEqual(testCase, actual, expected);
end

function testZeropad1DReal(testCase)
    N = 32;
    domain = PSDomain(setup1dX(N), false, false);
    
    f = domain.x{1};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    
    expected = zeros(48,1);
    expected(1:N) = f;
    
    verifyEqual(testCase, actual, expected);
end

function testZeropad2DReal(testCase)
    N = 32;
    domain = PSDomain(setup2dX(N), false, false);
    
    f = domain.x{1} + domain.x{2};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    expected = zeros(48);
    expected(1:N,1:N/2) = f(:,1:N/2);
    expected(1:N,ratio*N-N/2+1:ratio*N) = f(:,1+N/2:end);
    
    verifyEqual(testCase, actual, expected);
end

function testTruncation1DReal(testCase)
    N = 3*32;
    domain = PSDomain(setup1dX(N), false, false);
    
    f = domain.x{1};
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    
    expected = f(1:N*2/3);
    
    verifyEqual(testCase, actual, expected);
end

function testTruncation2DReal(testCase)
    N = 3*32;
    domain = PSDomain(setup2dX(N), false, false);
    
    f = domain.x{1} + domain.x{2};
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    expected = f(1:N*2/3, [1:N*1/3, N*2/3+1:N]);
    
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
    coarseDomain = PSDomain(setup1dX(N), true);
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

function testConvDealiasing2(testCase)
    N = 2^5; M = 2^6;
    coarseDomain = PSDomain(setup1dX(N), true);
    fineDomain = PSDomain(setup1dX(M), false);
    
    actual = eval(coarseDomain, N);
    expected = eval(fineDomain, N);
    
%     hold on
%     plot(abs(actual))
%     plot(abs(expected([1:end/4,3/4*end+1:end])))
%     plot(angle(actual))
%     plot(angle(expected([1:end/4,3/4*end+1:end])))
    
    verifyEqual(testCase, actual, expected([1:end/4,3/4*end+1:end]), 'AbsTol', 1e-13);
    
    function out = eval(domain, N)
        % Remove offset between the two resolutions by zeroing x
        u = 1 + sin(N/2 * pi * (domain.x{1} - domain.x{1}(1)));
        v = cos((N/2 + 2) * pi * (domain.x{1} - domain.x{1}(1)));
        uhat = domain.fft(u);
        vhat = domain.fft(v);
        out = domain.multiply(uhat, vhat);
    end
end

function testConvDealiasingReal(testCase)
    N = 2^5; M = 2^6;
    coarseDomain = PSDomain(setup1dX(N), true, false);
    fineDomain = PSDomain(setup1dX(M), false, false);
    
    actual = eval(coarseDomain, N);
    expected = eval(fineDomain, N);
    
%     hold on
%     plot(real(actual))
%     plot(real(expected(1:end/2)))
%     plot(imag(actual))
%     plot(imag(expected))
    
    verifyEqual(testCase, actual, expected(1:end/2), 'AbsTol', 1e-13);
    
    function out = eval(domain, N)
        u = 1 + sin(N/2 * pi * (domain.x{1} - domain.x{1}(1)));
        v = cos((N/2 + 2) * pi * (domain.x{1} - domain.x{1}(1)));
        uhat = domain.fft(u);
        vhat = domain.fft(v);
        out = domain.multiply(uhat, vhat);
    end
end

function testConvDealiasingHigherPowers(testCase)
    N = 2^7; M = 2^8;
    coarseDomain = PSDomain(setup1dX(N), true);
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
    coarseDomain = PSDomain(setup1dX(N), true);
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
