function out = diff_ps_2d(x, y, deg)
    %DIFF_PSEUDO_SPECTRAL Uses the pseudo-spectral method to differentiate
    %   Detailed explanation goes here
    suppression = 1e-13;
    
     % Transform into fourier space
    fy = fft2(y);
    % Prior suppression
    fy(abs(fy)<suppression) = 0;
    
    out = zeros([size(y),length(deg)]);
    N = cellfun(@(x) length(x)/2,x);
    L = cellfun(@(x) x(end),x);
    k = cell(2,1);
    for j = 1:length(x)
        k{j} = repmat([0:N(j)-1, 0, 1-N(j):-1]',1,2) * 2*pi./L(j);
    end
    
    for di = 1:length(deg)
        % Apply pseudo-spectral differentiation
        % Transform back into real space
        out(:,:,di,:) = ifftn(((1i*k{1}).^deg(1,di)*(1i*k{2}').^deg(2,di)).*fy);
    end
end
