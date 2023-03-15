%% Inverse Sympletic Finite Fourier Transform
% Author: Yunbo Hu
function X_out = ISFFT(X)
    m = size(X,1);
    n = size(X,2);
    X_out = ifft(fft(X,m).',n).' * sqrt(n/m); %ISFFT -> DFT in delay domain and IDFT in Doppler domain
end

