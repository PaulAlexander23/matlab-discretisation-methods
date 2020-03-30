function tests = testCTFDDomain()
    tests = functiontests(localfunctions);
end

function test1DDiff(testCase)
    N = 1024;
    x = {linspace(2*pi/N, 2*pi, N)'};
    degrees = [1, 2];
    accuracy = 4;
    domain = CTFDDomain(x, degrees, accuracy);
    
    f = exp(cos(domain.x{1}));

    expected = -sin(domain.x{1}) .* exp(cos(domain.x{1}));
    actual = domain.diff(f, 1);
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-5, 'AbsTol', 1e-7);

    expected = sin(domain.x{1}).^2 .* exp(cos(domain.x{1})) - cos(domain.x{1}) .* exp(cos(domain.x{1}));
    actual = domain.diff(f, 2);
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-5, 'AbsTol', 1e-7);
end

function test1DDiffDegree1Convergence(testCase)
    N = ceil(logspace(1, 3));
    errorVector = zeros(length(N), 1);
    
    degree = 1;
    for accuracy = 2:2:8
        for n = 1:length(N)
            x = {linspace(2*pi/N(n), 2*pi, N(n))'};
            domain = CTFDDomain(x, degree, accuracy);
            
            f = exp(cos(domain.x{1}));
            
            expected = -sin(domain.x{1}) .* exp(cos(domain.x{1}));
            actual = domain.diff(f, 1);
            
            errorVector(n) = max(abs(actual - expected));
        end

        viableIndices = logical((errorVector < 10^(-1.5)) .* (errorVector > 1e-11) .* (N' > 20));
        actualConvergence = -mean(gradient(log10(errorVector(viableIndices)), log10(N(viableIndices))));

        verifyEqual(testCase, actualConvergence, accuracy, 'RelTol', 3e-2);
    end
end

function test1DDiffDegree2Convergence(testCase)
    N = ceil(logspace(1, 3));
    errorVector = zeros(length(N), 1);
    
    degree = 2;
    for accuracy = 2:2:8
        for n = 1:length(N)
            x = {linspace(2*pi/N(n), 2*pi, N(n))'};
            domain = CTFDDomain(x, degree, accuracy);
            
            f = exp(cos(domain.x{1}));
            
            expected = sin(domain.x{1}).^2 .* exp(cos(domain.x{1})) - cos(domain.x{1}) .* exp(cos(domain.x{1}));
            actual = domain.diff(f, 2);
            
            errorVector(n) = max(abs(actual - expected));
        end

        viableIndices = logical((errorVector < 10^(-1.5)) .* (errorVector > 10^(-9.5)) .* (N' > 20));
        actualConvergence = -mean(gradient(log10(errorVector(viableIndices)), log10(N(viableIndices))));

        verifyEqual(testCase, actualConvergence, accuracy, 'RelTol', 3e-2);
    end
end

function test1DDiffMatInversion(testCase)
    N = 1024;
    x = {linspace(2*pi/N, 2*pi, N)'};
    degrees = [1, 2];
    accuracy = 4;
    domain = CTFDDomain(x, degrees, accuracy);
    [DA, DB] = domain.diffMat(1);
    [DA2, D2B] = domain.diffMat(2);
    
    f = exp(cos(domain.x{1}));

    expected = domain.diff(f, 1);
    actual = DA \ (DB * f);
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-5, 'AbsTol', 1e-7);

    expected = domain.diff(f, 2);
    actual = DA2 \ (D2B * f);
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-5, 'AbsTol', 1e-7);
end

function test1DDiffMat(testCase)
    N = 1024;
    x = {linspace(2*pi/N, 2*pi, N)'};
    degrees = [1, 2];
    accuracy = 4;
    domain = CTFDDomain(x, degrees, accuracy);
    D = domain.diffMat(1);
    D2 = domain.diffMat(2);
    
    f = exp(cos(domain.x{1}));

    expected = domain.diff(f, 1);
    actual = D * f;
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-5, 'AbsTol', 1e-7);

    expected = domain.diff(f, 2);
    actual = D2 * f;
    verifyEqual(testCase, actual, expected, 'RelTol', 1e-5, 'AbsTol', 1e-7);
end
