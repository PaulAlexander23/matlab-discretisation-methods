function d = finiteDifferenceCoefficients(M, N, x0, a)
    %Implementation of algorith found in Generation of Finite Difference
    %Fomulas on Arbitrarily Spaced Grids, by Bengt Fornberg, 1988
    d = zeros(M + 2, N + 1, N + 1);
    d(2,1,1) = 1;
    c1 = 1;
    
    for n = 1:N
        c2 = 1;
        for nu = 0:n - 1
            c3 = a(n+1) - a(nu+1);
            c2 = c2 * c3;
            if n <= M, d(n+2, n, nu+1) = 0; end
            for m = 0:min(n, M)
                d(m+2,n+1,nu+1) = ((a(n+1) - x0) * d(m+2,n,nu+1) - m * d(m+1,n,nu+1))/c3;
            end
        end
        for m = 0:min(n, M)
            d(m+2,n+1,n+1) = (c1/c2) * (m * d(m+1,n,n) - (a(n) - x0) * d(m+2,n,n));
        end
        c1 = c2;
    end
end