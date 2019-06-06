function v = dealiasZeroPad(u, p)
    %DEALIASZEROPAD Dealiasing via zero padding as in Fei Lu's report.
    %   u - surface to be dealiased
    %   p - nonlinearity ie (u^p)_x
    [N, tN] = size(u);
    K = (p + 1) * N / 2;
    upad = zeros(K,tN);
    indupad = [1:N/2, K-N/2+1:K];
    upad(indupad, :) = u;
    temp = fft(real(ifft(upad).^2));
    temp = K/N * temp;
    v = temp(indupad);
    v(N/2 + 1) = 0;
end