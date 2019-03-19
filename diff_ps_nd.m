function out = diff_ps_nd(x, y, deg)
    %DIFF_PSEUDO_SPECTRAL Uses the pseudo-spectral method to differentiate
    %   Detailed explanation goes here
    suppression = 1e-13;
    ndeg = length(deg);
    ndim = size(deg,1);
    out = zeros([size(y),ndeg]);
    N = size(fy,1)/2;
    % Determine k in matlab form
    %k = fftshift([0,-N+1:N-1])';
    k = repmat([0:N-1, 0, 1-N:-1]',1,ndim) * 2*pi./x(end,:);
    
    % Transform into fourier space
    fy = fftn(y);
    % Prior suppression
    fy(abs(fy)<suppression) = 0;
    
    matrixdims = repmat(':',1,ndim-1);
    for di = 1:ndeg
        % Apply pseudo-spectral differentiation
        % Transform back into real space
        iK = eye(size(y));
        for j = 1:ndim
            error('Not implemented')
            iK = iK * (1i*k(:,j)).^deg(di,j);% need to transpose each vector
        end
        out(matrixdims,di) = ifftn(iK.*fy);
    end
end
