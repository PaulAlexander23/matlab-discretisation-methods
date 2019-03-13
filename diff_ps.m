function out = diff_ps(x, y, degree)
    %DIFF_PSEUDO_SPECTRAL Uses the pseudo-spectral method to differentiate
    %   Detailed explanation goes here
    suppression = 1e-13;
    
    % Transform into fourier space
    out = fft(y);
    
    N = size(out,1)/2;
    
    % Determine k in matlab form
    %k = fftshift([0,-N+1:N-1])';
    k = [0:N-1, 0, 1-N:-1]' * 2*pi/x(end);
    
    % Prior suppression
    out(abs(out)<suppression) = 0;
    
    % Apply pseudo-spectral differentiation
    out = (1i*k).^degree.*out;
    
    % Posterior suppression
    % dyF(abs(dyF) < suppression*N*2) = 0 ;
    % dyF(abs(dyF) < suppression*max(abs(dyF))) = 0 ;
    
    % Transform back into real space
    out = ifft(out);
    
end
