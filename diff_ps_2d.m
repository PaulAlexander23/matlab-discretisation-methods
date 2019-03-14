function out = diff_ps_2d(x, y, deg)
    %DIFF_PSEUDO_SPECTRAL Uses the pseudo-spectral method to differentiate
    %   Detailed explanation goes here
    suppression = 1e-13;
    
     % Transform into fourier space
    fy = fft2(y);
    % Prior suppression
    fy(abs(fy)<suppression) = 0;
    
    out = zeros([size(y),length(deg)]);
    N = size(fy,1)/2;
    k = repmat([0:N-1, 0, 1-N:-1]',1,2) * 2*pi./x(end,:);
    
    for di = 1:length(deg)
        % Apply pseudo-spectral differentiation
        % Transform back into real space
        out(:,:,di) = ifftn(((1i*k(:,1)).^deg(di,1)*(1i*k(:,2)').^deg(di,2)).*fy);
    end
end
