function v = mydealiasTrunc(u, p)
    %DEALIASTRUNC Dealiasing via zero padding as in Fei Lu's report.
    %   u - surface to be dealiased
    %   p - nonlinearity ie (u^p)_x
    r = (p + 1)/2;
    upad = trunc(u,1/r);
    temp = fft(real(ifft(upad).^2));
    temp = r * temp;
    v = zeropad(temp,r);
    %v(end/2 + 1) = 0;
end