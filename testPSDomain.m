function tests = testPSDomain()
    tests = functiontests(localfunctions);
end

function testConstructor1DDefault(testCase)
    domain = PSDomain(setup1dX(2^8));
    
    verifyEqual(testCase, domain.length, 1);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        [0:2^7 - 1, 0, 1 - 2^7:-1]' * 2*pi);
    verifyEqual(testCase, domain.antialiasing, false);
    verifyEqual(testCase, domain.complex, true);
    verifyEqual(testCase, domain.suppression, eps);
end

function testConstructor1DAntiAliasingFlag(testCase)
    domain = PSDomain(setup1dX(2^8), false);
    verifyEqual(testCase, domain.antialiasing, false);

    domain = PSDomain(setup1dX(2^8), true);
    verifyEqual(testCase, domain.antialiasing, true);
end

function testConstructor1DComplexFlag(testCase)
    domain = PSDomain(setup1dX(2^8), true, false);
    verifyEqual(testCase, domain.complex, false);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        (0:2^7 - 1)' * 2*pi);

    domain = PSDomain(setup1dX(2^8), true, true);
    verifyEqual(testCase, domain.complex, true);

    domain = PSDomain(setup1dX(2^8), false, false);
    verifyEqual(testCase, domain.complex, false);
    verifyEqual(testCase, domain.wavenumber{1}, ...
        (0:2^7 - 1)' * 2*pi);

    domain = PSDomain(setup1dX(2^8), false, true);
    verifyEqual(testCase, domain.complex, true);
end

function testReshapeToVectorReal(testCase)
    N = 8;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N), linspace(L/N,L,N)}, false, false);
    y = domain.fft(cos(domain.x{1}) + cos(domain.x{2}));

    expected = reshape(y, [], 1);
    actual = domain.reshapeToVector(y);

    verifyEqual(testCase, actual, expected);

    z = [y, 2*y; 3*y, 4*y];
    yv = reshape(y, [], 1);

    expected = [yv, 2*yv; 3*yv, 4*yv];
    actual = domain.reshapeToVector(z);

    verifyEqual(testCase, actual, expected);
end

function testReshapeToDomainReal(testCase)
    N = 8;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N), linspace(L/N,L,N)}, false, false);
    Y = (1:N*(N/2))';

    expected = (1:N/2)' + (0:N-1)*N/2;
    actual = domain.reshapeToDomain(Y);

    verifyEqual(testCase, actual, expected);

    Y = (1:N*(N/2))';
    y = (1:N/2)' + (0:N-1)*N/2;

    expected = [y, 2*y; 3*y, 4*y];
    actual = domain.reshapeToDomain([Y, 2*Y; 3*Y, 4*Y]);

    verifyEqual(testCase, actual, expected);
end

function testFFT1D(testCase)
    N = 2^7;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)});

    for m = 1:N/2-1
        y = cos(m * 2*pi*(domain.x{1} - domain.x{1}(1)));
        
        actual = domain.fft(y);
        expected = zeros(N, 1);
        expected(m+1) = 1;
        expected(N-m+1) = 1;
        
        verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
    end
end

function testFFT1DReal(testCase)
    N = 2^7;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)}, false, false);

    for m = 1:N/2-1
        y = cos(m * 2*pi*(domain.x{1} - domain.x{1}(1)));
        
        actual = domain.fft(y);
        expected = zeros(N/2, 1);
        expected(m+1) = 1;
        
        verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
    end
end

function testFFTMultipleSurfaces(testCase)
    N = 128;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N),linspace(L/N,L,N)}, false, true);
    y = cos(domain.x{1} + domain.x{2});
    z = [y, 2*y; 3*y, 4*y];

    fy = domain.fft(y);
    expected = [fy, 2*fy; 3*fy, 4*fy];
    actual = domain.fft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 2*eps);
end

function testFFTMultipleSurfacesReal(testCase)
    N = 128;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N),linspace(L/N,L,N)}, false, false);
    y = cos(domain.x{1} + domain.x{2});
    z = [y, 2*y; 3*y, 4*y];

    fy = domain.fft(y);
    expected = [fy, 2*fy; 3*fy, 4*fy];
    actual = domain.fft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 2*eps);
end

function testIFFT1D(testCase)
    N = 2^4;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)});
    
    for m = 1:N/2-1
        y = zeros(N, 1);
        y(m+1) = 1;
        y(N-m+1) = 1;
        
        expected = cos(m * 2*pi*(domain.x{1} - domain.x{1}(1)));
        actual = domain.ifft(y);
    
        verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14);
    end
end

function testIFFT1DReal(testCase)
    N = 2^4;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)}, false, false);
    
    for m = 1:N/2-1
        y = zeros(N/2, 1);
        y(m+1) = 1;
        
        expected = cos(m * 2*pi*(domain.x{1} - domain.x{1}(1)));
        actual = domain.ifft(y);
    
        verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14);
    end
end

function testIFFT1DMultipleSurfaces(testCase)
    N = 2^4;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)});
    
    y = zeros(N, 1);
    y(2) = 1;
    y(N) = 1;
    
    z = [y, 2*y; 3*y, 4*y];
    Y = cos(2*pi*(domain.x{1} - domain.x{1}(1)));

    expected = [Y, 2*Y; 3*Y, 4*Y];
    actual = domain.ifft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14);
end

function testIFFT1DMultipleSurfacesReal(testCase)
    N = 2^4;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)}, false, false);
    
    y = zeros(N/2, 1);
    y(2) = 1;

    z = [y, 2*y; 3*y, 4*y];
    Y = cos(2*pi*(domain.x{1} - domain.x{1}(1)));

    expected = [Y, 2*Y; 3*Y, 4*Y];
    actual = domain.ifft(z);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-14, 'AbsTol', 1e-14);
end

function testDiff1D(testCase)
    N = 2^7;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N)});

    for m = 1:N/2-1
        y = domain.fft(cos(m * domain.x{1}));
        
        actual = domain.diff(y, 1);
        expected = domain.fft(-m * sin(m * domain.x{1}));
        
        verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);

        actual = domain.diff(y, 2);
        expected = domain.fft(-m^2 * cos(m * domain.x{1}));
        
        verifyEqual(testCase, actual, expected, 'RelTol', 1e-10, 'AbsTol', 1e-10);
    end
end

function testDiff1DReal(testCase)
    N = 2^7;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N)}, false, false);

    for m = 1:N/2-1
        y = domain.fft(cos(m * domain.x{1}));
        
        actual = domain.diff(y, 1);
        expected = domain.fft(-m * sin(m * domain.x{1}));
        
        verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);

        actual = domain.diff(y, 2);
        expected = domain.fft(-m^2 * cos(m * domain.x{1}));
        
        verifyEqual(testCase, actual, expected, 'RelTol', 1e-10, 'AbsTol', 1e-10);
    end
end

function testDiff2D(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)});

    for m = 1:M/2-1
        for n = 1:N/2-1
            y = domain.fft(cos(m * domain.x{1} + n * domain.x{2}));
            
            actual = domain.diff(y, [1, 0]');
            expected = domain.fft(-m * sin(m * domain.x{1} + n * domain.x{2}));
            
            verifyEqual(testCase, actual, expected, 'RelTol', 1e-11, 'AbsTol', 1e-11);

            actual = domain.diff(y, [0, 1]');
            expected = domain.fft(-n * sin(m * domain.x{1} + n * domain.x{2}));
            
            verifyEqual(testCase, actual, expected, 'RelTol', 1e-11, 'AbsTol', 1e-11);
        end
    end
end

function testDiff2DReal(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)}, false, false);

    for m = 1:M/2-1
        for n = 1:N/2-1
            y = domain.fft(cos(m * domain.x{1} + n * domain.x{2}));
            
            actual = domain.diff(y, [1, 0]');
            expected = domain.fft(-m * sin(m * domain.x{1} + n * domain.x{2}));
            
            verifyEqual(testCase, actual, expected, 'RelTol', 1e-11, 'AbsTol', 1e-11);

            actual = domain.diff(y, [0, 1]');
            expected = domain.fft(-n * sin(m * domain.x{1} + n * domain.x{2}));
            
            verifyEqual(testCase, actual, expected, 'RelTol', 1e-11, 'AbsTol', 1e-11);
        end
    end
end

function testDiff2DMultipleSurfaces(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)});

    y = domain.fft(cos(domain.x{1} + domain.x{2}));
    dy = domain.fft(-sin(domain.x{1} + domain.x{2}));
    
    actual = domain.diff([y, 2*y; 3*y, 4*y], [1, 0]');
    expected = [dy, 2*dy; 3*dy, 4*dy];
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-11, 'AbsTol', 1e-11);

    actual = domain.diff([y, 2*y; 3*y, 4*y], [0, 1]');
    expected = [dy, 2*dy; 3*dy, 4*dy];
    
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-11, 'AbsTol', 1e-11);
end

function testGetDiffMatrix1D(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)});
    
    expected = spdiags(1i * [0:N/2-1,0,-N/2+1:-1]', 0, N, N);
    actual = domain.diffMat(1);
    
    verifyEqual(testCase, actual, expected);

    expected = spdiags(- [0:N/2-1,0,-N/2+1:-1]' .^ 2, 0, N, N);
    actual = domain.diffMat(2);
    
    verifyEqual(testCase, actual, expected);
end

function testGetDiffMatrix2D(testCase)
    M = 2^4;
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)});
    
    expected = kron(speye(N), spdiags(1i * [0:M/2-1,0,-M/2+1:-1]', 0, M, M));
    actual = domain.diffMat([1, 0]');
    
    verifyEqual(testCase, actual, expected);

    expected = kron(spdiags(1i * [0:N/2-1,0,-N/2+1:-1]', 0, N, N), speye(M));
    actual = domain.diffMat([0, 1]');
    
    verifyEqual(testCase, actual, expected);
end

function testMultiply1D(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)},false);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            expected = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1})), domain.fft(cos(n * domain.x{1})));

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testMultiply1DDealiased(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)}, true);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            expected = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1})), domain.fft(cos(n * domain.x{1})));

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end

    for n = 1:N/2 - 1
        for m = N/2 - n:N/2 - 1 
            expected = domain.fft(cos((m - n) * domain.x{1}));
            actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1})), domain.fft(cos(n * domain.x{1})));

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testZeropad1D(testCase)
    N = 32;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)});
    
    f = domain.x{1};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    
    expected = zeros(N * ratio, 1);
    expected(1:N/2) = f(1:N/2);
    expected(ratio*N-N/2+1:ratio*N) = f(1+N/2:end);
    
    verifyEqual(testCase, actual, expected);
end

function testZeropad2D(testCase)
    M = 64;
    N = 32;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)});
    
    f = domain.x{1} + domain.x{2};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    
    expected = zeros([M, N] * ratio);
    expected(1:M/2,1:N/2) = f(1:M/2,1:N/2);
    expected(1:M/2,ratio*N-N/2+1:ratio*N) = f(1:M/2,1+N/2:end);
    expected(ratio*M-M/2+1:ratio*M,1:N/2) = f(1+M/2:end,1:N/2);
    expected(ratio*M-M/2+1:ratio*M,ratio*N-N/2+1:ratio*N) = f(1+M/2:end,1+N/2:end);
    
    verifyEqual(testCase, actual, expected);
end

function testZeropad1DReal(testCase)
    N = 32;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)}, false, false);
    
    f = domain.x{1};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    
    expected = zeros(N * ratio,1);
    expected(1:N) = f;
    
    verifyEqual(testCase, actual, expected);
end

function testZeropad2DReal(testCase)
    M = 64;
    N = 32;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)}, false, false);
    
    f = domain.x{1} + domain.x{2};
    ratio = 3/2;
    
    actual = domain.zeropad(f, ratio);
    expected = zeros([M, N] * ratio);
    expected(1:M,1:N/2) = f(1:M,1:N/2);
    expected(1:M,ratio*N-N/2+1:ratio*N) = f(1:M,1+N/2:end);

    verifyEqual(testCase, actual, expected);
end

function testTrunc1D(testCase)
    N = 3*32;
    domain = PSDomain(setup1dX(N));
    
    f = [0:N/2-1,0,-N/2+1:-1]';
    ratio = 2/3;
    
    expected = [0:ratio*N/2-1,0,-ratio*N/2+1:-1]';
    actual = domain.trunc(f, ratio);
    
    verifyEqual(testCase, actual, expected);
end

function testTrunc2D(testCase)
    N = 3*32;
    domain = PSDomain(setup2dX(N));
    
    f = repmat([0:N/2-1,0,-N/2+1:-1]', 1, N);
    ratio = 2/3;
    
    expected = repmat([0:ratio*N/2-1,0,-ratio*N/2+1:-1]', 1, N*ratio);
    actual = domain.trunc(f, ratio);
    
    verifyEqual(testCase, actual, expected);
end

function testTrunc1DReal(testCase)
    N = 3*32;
    domain = PSDomain(setup1dX(N), false, false);
    
    f = (0:N-1)';
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    
    expected = (0:N*ratio-1)';
    
    verifyEqual(testCase, actual, expected);
end

function testTrunc2DReal(testCase)
    N = 3*32;
    domain = PSDomain(setup2dX(N), false, false);
    
    f = repmat((0:N-1)', 1, N);
    ratio = 2/3;
    
    actual = domain.trunc(f, ratio);
    expected = repmat((0:N*ratio-1)', 1, N*ratio);
    
    verifyEqual(testCase, actual, expected);
end

function testScaling1D(testCase)
    n = 2^8;
    domain = PSDomain(setup1dX(n));
    amplitude = 2;
    y = amplitude * cos(2 * pi * domain.x{1});
    fy = domain.fft(y);
    verifyEqual(testCase, max(real(fy)), amplitude, 'RelTol', 1e-3);
end

function testScaling2D(testCase)
    n = 2^8;
    domain = PSDomain(setup2dX(n));
    amplitude = 2;
    y = amplitude * cos(2 * pi * (domain.x{1} + domain.x{2}));
    fy = domain.fft(y);
    verifyEqual(testCase, max(max(real(fy))), amplitude, 'RelTol', 1.3e-3);
end

%% Functions
function x = setup1dX(n)
    x = {linspace(1/n, 1, n)'};
end

function x = setup2dX(m,n)
    if nargin < 2
        n = m;
    end
    x = {linspace(1/m,1,m)';linspace(1/n,1,n)};
end
