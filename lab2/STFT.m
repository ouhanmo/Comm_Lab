function [S, f, t] = STFT(x, window, Noverlap, Nfft, fs)
f = fs/Nfft:fs/Nfft:fs;
f = [f(floor(1/2*(Nfft)+1):end)-fs,f(1:floor(1/2*(Nfft)))];
step = length(window)-Noverlap;
t = 1/fs*(1:step:length(x));
S = zeros(Nfft, length(t));

x_tap = [zeros(1,floor(length(window)/2)) x zeros(1,floor(length(window)/2)+1)];

for seg = 0:length(t)-1
    segment = x_tap(1+seg*step:seg*step+length(window));
    circshift(fft(segment.*window,Nfft),floor(1/2*(Nfft+1)));
    S(1:end,seg+1) = circshift(fft(segment.*window,Nfft),floor(1/2*(Nfft)));
end

