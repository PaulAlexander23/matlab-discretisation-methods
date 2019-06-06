function v = zeropad(u, r)
    %ZEROPAD Zero pad vector or matrix in the first dimension
    %   u - vector or matrix
    %   r - extention ratio
    [N, M] = size(u);
    v = zeros(N * r, M);
    ind = [1:N/2, N * (r - 0.5) + 1:N * r];
    v(ind,:) = u;
end