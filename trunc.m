function v = trunc(u, r)
    %TRUNC Truncate vector or matrix in the first dimension
    %   u - vector or matrix
    %   r - truncation ratio
    N = size(u,1);
    ind = [1:r * N/2, N * (1 - r/2) + 1:N];
    v = u(ind,:);
end