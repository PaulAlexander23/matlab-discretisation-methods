function tests = testPSDomain()
    tests = functiontests(localfunctions);
end

function testConstructorDefault(testCase)
    domain = PSDomain(setup1dX(2^8));

    verifyEqual(testCase, domain.antialiasing, false);
    verifyEqual(testCase, domain.complex, true);
    verifyEqual(testCase, domain.suppression, 1e-15);

    domain = PSDomain(setup2dX(2^8));

    verifyEqual(testCase, domain.antialiasing, false);
    verifyEqual(testCase, domain.complex, true);
    verifyEqual(testCase, domain.suppression, 1e-15);
end

function testLength(testCase)
    N = 2^5;
    L = 2*pi;

    domain = PSDomain({linspace(L/N, L, N)});
    verifyEqual(testCase, domain.length, L);

    domain = PSDomain({linspace(L/N/2, L - L/N/2, N)});
    verifyEqual(testCase, domain.length, L, 'RelTol', eps);

    domain = PSDomain({linspace(L/N, L, N), linspace(L/N, L, N)});
    verifyEqual(testCase, domain.length, [L; L]);

    domain = PSDomain({linspace(L/N/2, L - L/N/2, N), linspace(L/N/2, L - L/N/2, N)});
    verifyEqual(testCase, domain.length, [L; L], 'RelTol', eps);
end

function testWavenumber(testCase)
    N = 2^5;
    L = 2*pi;

    domain = PSDomain({linspace(L/N, L, N)});
    verifyEqual(testCase, domain.wavenumber, {[0:N/2-1, 0, -N/2+1:-1]'});

    domain = PSDomain({linspace(L/N, L, N)}, false, false);
    verifyEqual(testCase, domain.wavenumber, {(0:N/2-1)'});

    domain = PSDomain({linspace(L/N, L, N), linspace(L/N, L, N)});
    verifyEqual(testCase, domain.wavenumber, {[0:N/2-1, 0, -N/2+1:-1]', [0:N/2-1, 0, -N/2+1:-1]'});

    domain = PSDomain({linspace(L/N, L, N), linspace(L/N, L, N)}, false, false);
    verifyEqual(testCase, domain.wavenumber, {(0:N/2-1)', [0:N/2-1, 0, -N/2+1:-1]'});
end

function testSuppression(testCase)
    N = 2^5;
    L = 2*pi;

    domain = PSDomain({linspace(L/N, L, N)});
    verifyEqual(testCase, domain.suppression, 1e-15);

    domain = PSDomain({linspace(L/N, L, N)}, false, false, 1e-5);
    verifyEqual(testCase, domain.suppression, 1e-5);
end

function testAntiAliasingFlag(testCase)
    domain = PSDomain(setup1dX(2^8), false);
    verifyEqual(testCase, domain.antialiasing, false);

    domain = PSDomain(setup1dX(2^8), true);
    verifyEqual(testCase, domain.antialiasing, true);
end

function testComplexFlag(testCase)
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

function testScaling1D(testCase)
    n = 2^8;
    domain = PSDomain(setup1dX(n));
    amplitude = 2;

    y = amplitude * ones(n,1);
    fy = domain.fft(y);

    verifyEqual(testCase, max(real(fy)), amplitude, 'RelTol', eps);

    y = amplitude * cos(2 * pi * (domain.x{1} - domain.x{1}(1)));
    fy = domain.fft(y);

    verifyEqual(testCase, max(real(fy)), amplitude, 'RelTol', eps);

    y = amplitude * cos(2 * pi * (n/2)  * (domain.x{1} - domain.x{1}(1)));
    fy = domain.fft(y);

    verifyEqual(testCase, max(real(fy)), amplitude, 'RelTol', eps);
end

function testScaling2D(testCase)
    n = 2^8;
    domain = PSDomain(setup2dX(n));
    amplitude = 2;

    y = amplitude * ones(n);
    fy = domain.fft(y);
    verifyEqual(testCase, max(max(real(fy))), amplitude, 'RelTol', eps);

    y = amplitude * cos(2 * pi * (domain.x{1} - domain.x{1}(1) + domain.x{2} - domain.x{2}(1)));
    fy = domain.fft(y);
    verifyEqual(testCase, max(max(real(fy))), amplitude, 'RelTol', eps);
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

    for m = 0:N/2-1
        y = cos(m * 2*pi*(domain.x{1} - domain.x{1}(1)));

        actual = domain.fft(y);
        expected = zeros(N, 1);
        expected(m+1) = 1;
        if m>0
            expected(N-m+1) = 1;
        end

        verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
    end
end

function testFFT1DReal(testCase)
    N = 2^7;
    L = 1;
    domain = PSDomain({linspace(L/N,L,N)}, false, false);

    for m = 0:N/2-1
        y = cos(m * 2*pi*(domain.x{1} - domain.x{1}(1)));

        actual = domain.fft(y);
        expected = zeros(N/2, 1);
        expected(m+1) = 1;

        verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
    end
end

function testFFT1DMultipleSurfaces(testCase)
    N = 128;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N)}, false, true);
    y = cos(domain.x{1});
    z = [y, 2*y; 3*y, 4*y];

    fy = domain.fft(y);
    expected = [fy, 2*fy; 3*fy, 4*fy];
    actual = domain.fft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 2*eps);
end

function testFFT1DMultipleSurfacesReal(testCase)
    N = 128;
    L = 2*pi;
    domain = PSDomain({linspace(L/N,L,N)}, false, false);
    y = cos(domain.x{1});
    z = [y, 2*y; 3*y, 4*y];

    fy = domain.fft(y);
    expected = [fy, 2*fy; 3*fy, 4*fy];
    actual = domain.fft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 2*eps);
end

function testFFT2D(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M,L,M),linspace(L/N,L,N)});

    for m = 1:M/2-1
        for n = 1:M/2-1
            y = cos(m * (domain.x{1} - domain.x{1}(1)) + n * (domain.x{2} - domain.x{2}(1)));

            actual = domain.fft(y);
            expected = zeros(M, N);
            expected(m+1, n+1) = 1;
            expected(M-m+1, N-n+1) = 1;

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testFFT2DReal(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M,L,M),linspace(L/N,L,N)}, false, false);

    for m = 1:M/2-1
        for n = 1:M/2-1
            y = cos(m * (domain.x{1} - domain.x{1}(1)) + n * (domain.x{2} - domain.x{2}(1)));

            actual = domain.fft(y);
            expected = zeros(M/2, N);
            expected(m+1, n+1) = 1;

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testFFT2DMultipleSurfaces(testCase)
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

function testFFT2DMultipleSurfacesReal(testCase)
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

    for m = 0:N/2-1
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

function testIFFT2D(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M,L,M), linspace(L/N,L,N)});

    for m = 1:M/2-1
        for n = 1:N/2-1
            y = zeros(M, N);
            y(m+1, n+1) = 1;
            y(M-m+1, N-n+1) = 1;

            expected = cos(m * (domain.x{1} - domain.x{1}(1)) + n * (domain.x{2} - domain.x{2}(1)));
            actual = domain.ifft(y);

            verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        end
    end
end

function testIFFT2DReal(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M,L,M), linspace(L/N,L,N)}, false, false);

    for m = 1:M/2-1
        for n = 1:N/2-1
            y = zeros(M/2, N);
            y(m+1, n+1) = 1;

            expected = cos(m * (domain.x{1} - domain.x{1}(1)) + n * (domain.x{2} - domain.x{2}(1)));
            actual = domain.ifft(y);

            verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        end
    end
end

function testIFFT2DMultipleSurfaces(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M,L,M), linspace(L/N,L,N)});

    y = zeros(M, N);
    y(2,2) = 1;
    y(M,N) = 1;

    z = [y, 2*y; 3*y, 4*y];
    Y = cos((domain.x{1} - domain.x{1}(1)) + (domain.x{2} - domain.x{2}(1)));

    expected = [Y, 2*Y; 3*Y, 4*Y];
    actual = domain.ifft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
end

function testIFFT2DMultipleSurfacesReal(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M,L,M), linspace(L/N,L,N)}, false, false);

    y = zeros(M/2, N);
    y(2,2) = 1;

    z = [y, 2*y; 3*y, 4*y];
    Y = cos((domain.x{1} - domain.x{1}(1)) + (domain.x{2} - domain.x{2}(1)));

    expected = [Y, 2*Y; 3*Y, 4*Y];
    actual = domain.ifft(z);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
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

function testDiff2DMultipleSurfacesReal(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)}, false, false);

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

function testGetDiffMatrix1DReal(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)}, false, false);

    expected = spdiags(1i * (0:N/2-1)', 0, N/2, N/2);
    actual = domain.diffMat(1);

    verifyEqual(testCase, actual, expected);

    expected = spdiags(- (0:N/2-1)' .^ 2, 0, N/2, N/2);
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

function testGetDiffMatrix2DReal(testCase)
    M = 2^4;
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)}, false, false);

    expected = kron(speye(N), spdiags(1i * (0:M/2-1)', 0, M/2, M/2));
    actual = domain.diffMat([1, 0]');

    verifyEqual(testCase, actual, expected);

    expected = kron(spdiags(1i * [0:N/2-1,0,-N/2+1:-1]', 0, N, N), speye(M/2));
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

function testMultiply1DReal(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)},false,false);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            expected = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1})), domain.fft(cos(n * domain.x{1})));

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testMultiply1DDealiasedReal(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)}, true, false);

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

function testMultiply1DMultipleSurfaces(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)},false);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            f = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            g1 = domain.fft(cos(m * domain.x{1}));
            g2 = domain.fft(cos(n * domain.x{1}));

            expected = [f, 2*f; 3*f, 4*f];
            actual = 2 * domain.multiply([g1, 2*g1; 3*g1, 4*g1],[g2, g2; g2, g2]);

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testMultiply1DMultipleSurfacesDealiased(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)},true);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            f = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            g1 = domain.fft(cos(m * domain.x{1}));
            g2 = domain.fft(cos(n * domain.x{1}));

            expected = [f, 2*f; 3*f, 4*f];
            actual = 2 * domain.multiply([g1, 2*g1; 3*g1, 4*g1],[g2, g2; g2, g2]);

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testMultiply1DMultipleSurfacesReal(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)},false,false);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            f = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            g1 = domain.fft(cos(m * domain.x{1}));
            g2 = domain.fft(cos(n * domain.x{1}));

            expected = [f, 2*f; 3*f, 4*f];
            actual = 2 * domain.multiply([g1, 2*g1; 3*g1, 4*g1],[g2, g2; g2, g2]);

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testMultiply1DMultipleSurfacesDealiasedReal(testCase)
    N = 2^5;
    L = 2*pi;
    domain = PSDomain({linspace(L/N, L, N)},true,false);

    for n = 1:N/2 - 2
        for m = 1:N/2 - 1 - n
            f = domain.fft(cos((m - n) * domain.x{1}) + cos((m + n) * domain.x{1}));
            g1 = domain.fft(cos(m * domain.x{1}));
            g2 = domain.fft(cos(n * domain.x{1}));

            expected = [f, 2*f; 3*f, 4*f];
            actual = 2 * domain.multiply([g1, 2*g1; 3*g1, 4*g1],[g2, g2; g2, g2]);

            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
        end
    end
end

function testMultiply2D(testCase)
    M = 2^5;
    N = 2^4;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)}, false);

    for n = 1:M/2 - 2
        for m = 1:M/2 - 1 - n
            for b = 1:N/2 - 2
                for a = 1:N/2 - 1 - b
                    expected = domain.fft(cos((m - n) * domain.x{1} + (a - b) * domain.x{2}) + ...
                        cos((m + n) * domain.x{1} + (a + b) * domain.x{2}));
                    actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1} + a * domain.x{2})), domain.fft(cos(n * domain.x{1} + b * domain.x{2})));

                    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
                end
            end
        end
    end
end

function testMultiply2DDealiased(testCase)
    M = 2^4;
    N = 2^3;
    L = 2*pi;
    domain = PSDomain({linspace(L/M, L, M), linspace(L/N, L, N)}, true);

    for n = 1:M/2 - 2
        for m = 1:M/2 - 1 - n
            for b = 1:N/2 - 2
                for a = 1:N/2 - 1 - b
                    expected = domain.fft(cos((m - n) * domain.x{1} + (a - b) * domain.x{2}) + ...
                        cos((m + n) * domain.x{1} + (a + b) * domain.x{2}));
                    actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1} + a * domain.x{2})), domain.fft(cos(n * domain.x{1} + b * domain.x{2})));

                    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
                end
            end
        end
    end

    for n = 1:M/2 - 1
        for m = (M/2 - n):(M/2 - 1)
            for b = 1:N/2 - 1
                for a = N/2 - b:N/2 - 1
                    expected = domain.fft(cos((m - n) * domain.x{1} + (a - b) * domain.x{2}));
                    actual = 2 * domain.multiply(domain.fft(cos(m * domain.x{1} + a * domain.x{2})), domain.fft(cos(n * domain.x{1} + b * domain.x{2})));

                    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
                end
            end
        end
    end
end

function testZeroSmallModes(testCase)
    tolerance = 1e-6;

    domain = PSDomain({1:16});

    f = 10.^(-domain.x{1});

    expected = f .* (f >= tolerance);
    actual = domain.zeroSmallModes(f, tolerance);

    verifyEqual(testCase, actual, expected);

    domain.suppression = tolerance;

    actual = domain.zeroSmallModes(f);

    verifyEqual(testCase, actual, expected);
end

function testFilterOutShortWaves1D(testCase)
    domain = PSDomain({1:16});
    f = domain.x{1};

    expected = f;
    actual = domain.filterOutShortWaves(f,1);

    verifyEqual(testCase, actual, expected);

    expected = 0 * f;
    actual = domain.filterOutShortWaves(f,0);

    verifyEqual(testCase, actual, expected);

    expected = [1:5,0,0,0,0,0,0,12:16]';
    actual = domain.filterOutShortWaves(f,2/3);

    verifyEqual(testCase, actual, expected);
end

function testFilterOutShortWaves2D(testCase)
    domain = PSDomain({1:120,1:120});
    f = domain.x{1} + 120*(domain.x{2} - 1);

    expected = f;
    actual = domain.filterOutShortWaves(f,2);

    verifyEqual(testCase, actual, expected);

    expected = 0 * f;
    actual = domain.filterOutShortWaves(f,0);

    verifyEqual(testCase, actual, expected);

    newf = domain.filterOutShortWaves(f,2/3,"rectangle");
    expected = 6241;
    actual = numel(find(newf));

    verifyEqual(testCase, actual, expected);

    newf = domain.filterOutShortWaves(f,2/3,"ellipse");
    expected = 5013;
    actual = numel(find(newf));

    verifyEqual(testCase, actual, expected);

    newf = domain.filterOutShortWaves(f,2/3,"triangle");
    expected = 3121;
    actual = numel(find(newf));

    verifyEqual(testCase, actual, expected);
end

function testFilterOutShortWaves2DReal(testCase)
    domain = PSDomain({1:120,1:120}, false, false);
    f = domain.x{1}(1:end/2) + 120*(domain.x{2} - 1);

    expected = f;
    actual = domain.filterOutShortWaves(f,2);

    verifyEqual(testCase, actual, expected);

    expected = 0 * f;
    actual = domain.filterOutShortWaves(f,0);

    verifyEqual(testCase, actual, expected);

    newf = domain.filterOutShortWaves(f,2/3,"rectangle");
    expected = 3160;
    actual = numel(find(newf));

    verifyEqual(testCase, actual, expected);

    newf = domain.filterOutShortWaves(f,2/3,"ellipse");
    expected = 2546;
    actual = numel(find(newf));

    verifyEqual(testCase, actual, expected);

    newf = domain.filterOutShortWaves(f,2/3,"triangle");
    expected = 1600;
    actual = numel(find(newf));

    verifyEqual(testCase, actual, expected);
end

function testZeropad1D(testCase)
    domain = PSDomain({1:16});
    f = domain.x{1};

    expected = [f(1:8); zeros(16,1); f(9:16)];
    actual = domain.zeropad(f,2);

    verifyEqual(testCase, actual, expected);

    expected = [f(1:8); zeros(8,1); f(9:16)];
    actual = domain.zeropad(f,3/2);

    verifyEqual(testCase, actual, expected);
end

function testZeropad1DReal(testCase)
    domain = PSDomain({1:16},true,false);
    f = domain.x{1}(1:8);

    expected = [f; zeros(8,1)];
    actual = domain.zeropad(f,2);

    verifyEqual(testCase, actual, expected);

    expected = [f; zeros(4,1)];
    actual = domain.zeropad(f,3/2);

    verifyEqual(testCase, actual, expected);
end

function testZeropad2D(testCase)
    domain = PSDomain({1:16,1:32});
    f = domain.x{1} + domain.x{2};

    expectedSize = [16,32]*2;
    actual = domain.zeropad(f,2);

    verifySize(testCase, actual, expectedSize);

    expectedSize = [16,32]*3/2;
    actual = domain.zeropad(f,3/2);

    verifySize(testCase, actual, expectedSize);
end

function testZeropad2DReal(testCase)
    domain = PSDomain({1:16,1:32},true,false);
    f = domain.x{1}(1:end/2) + domain.x{2}(1:end/2);

    expectedSize = [8,16]*2;
    actual = domain.zeropad(f,2);

    verifySize(testCase, actual, expectedSize);

    expectedSize = [8,16]*3/2;
    actual = domain.zeropad(f,3/2);

    verifySize(testCase, actual, expectedSize);
end

function testTrunc1D(testCase)
    domain = PSDomain({1:16});
    f = domain.x{1};

    expected = [f(1:4); 0; f(14:16)];
    actual = domain.trunc(f,1/2);

    verifyEqual(testCase, actual, expected);

    expected = [f(1:6); 0; f(12:16)];
    actual = domain.trunc(f,3/4);

    verifyEqual(testCase, actual, expected);
end

function testTrunc1DReal(testCase)
    domain = PSDomain({1:16},true,false);
    f = domain.x{1}(1:8);

    expected = f(1:4);
    actual = domain.trunc(f,1/2);

    verifyEqual(testCase, actual, expected);

    expected = f(1:6);
    actual = domain.trunc(f,3/4);

    verifyEqual(testCase, actual, expected);
end

function testTrunc2D(testCase)
    domain = PSDomain({1:16,1:32});
    f = domain.x{1} + domain.x{2};

    expectedSize = [16,32]*2;
    actual = domain.trunc(f,2);

    verifySize(testCase, actual, expectedSize);

    expectedSize = [16,32]*3/2;
    actual = domain.trunc(f,3/2);

    verifySize(testCase, actual, expectedSize);
end

function testTrunc2DReal(testCase)
    domain = PSDomain({1:16,1:32},true,false);
    f = domain.x{1}(1:end/2) + domain.x{2}(1:end/2);

    expectedSize = [8,16]/2;
    actual = domain.trunc(f,1/2);

    verifySize(testCase, actual, expectedSize);

    expectedSize = [8,16]*3/4;
    actual = domain.trunc(f,3/4);

    verifySize(testCase, actual, expectedSize);
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
