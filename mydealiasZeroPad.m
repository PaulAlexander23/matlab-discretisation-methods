function v = mydealiasZeroPad(u, p)
    %DEALIASZEROPAD Dealiasing via zero padding as in Fei Lu's report.
    %   u - surface to be dealiased
    %   p - nonlinearity ie (u^p)_x
    r = (p + 1)/2;
    upad = zeropad(u,r);
    temp = fft(real(ifft(upad).^2));
    temp = r * temp;
    v = trunc(temp,1/r);
    %v(end/2 + 1) = 0;
end