%% Sympletic Finite Fourier Transform
% Author: Yunbo Hu
function X_out = SFFT(X)
    m = size(X,1);
    n = size(X,2);
    X_out = fft(ifft(X,m).',n).' * sqrt(m/n); %SFFT -> IDFT in delay domain and DFT in Doppler domain
end

