function output=ft(in)
% Performs fftshift(fft(fftshift(input)))
% 1-D forward FT
output = fftshift(fft(fftshift(in)));
