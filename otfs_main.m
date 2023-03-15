%% OTFS Sensing Simulation
% Author: Yunbo HU (Email: EdenHu1111@foxmail.com)
clc;clear;
%% Basic Paramaters
M = 1024;   %Number of data in delay domain
N = 120;   %Number of data in Doppler domain

c = 3e+8;
fc = 30e+9;
deltaf = 240e+3; % subcarrier interval 240kHz
T = 1/(deltaf);% period of pulse shaping filter
lambda = c/fc;

QAM_mod = 4;
bitPerSymbol = log2(QAM_mod);

eng_sqrt = (QAM_mod==2)+(QAM_mod~=2)*sqrt((QAM_mod-1)/6*(2^2)); % Average Power of a QAM symbol
SNR_dB = -10:5:20; 
NFrame = 100;

sigmaSquare = eng_sqrt^2 * exp(-SNR_dB*log(10)/10);

Nsymbolperframe = M*N; % Number of QAM symbols in a OTFS frame
Nbitsperframe = Nsymbolperframe * bitPerSymbol; % Number of bits in a OTFS frame
%% Infromation about the Target
% targetRange = [35,40];%target range 35m
targetRange = 35;
% targetSpeed = [20,30];%target speed 
targetSpeed = 30;%target speed 
delay = 2*targetRange/c;        kp = round(delay*M*deltaf);

doppler = 2*targetSpeed/lambda; lp = (doppler*N*T);

% %% basic characters
% P = eye(M*N);
% P = [P(M*N,:);P(1:M*N-1,:)]; %% shift matrix
% 
dd = 0:1:M*N-1;
d = exp(1j*2*pi/(M*N)*dd.');
% D = diag(d); 
% 
% H = zeros(M*N,M*N);
% 
% %% Channel Generation
% for iTarget = 1:length(targetRange)
%    H = H + P^kp(iTarget)*D^lp(iTarget)*(randn + 1j*randn)/sqrt(2);
% end
errorRange_TF = zeros(size(SNR_dB));
errorVelo_TF = zeros(size(SNR_dB));
errorRange_MF = zeros(size(SNR_dB));
errorVelo_MF = zeros(size(SNR_dB));
for ii = 1:length(SNR_dB)
    fprintf('SNR = %d\n',SNR_dB(ii));
    for jj = 1:NFrame
        data_info_bit = randi([0,1],Nbitsperframe,1);
        data_temp = bi2de(reshape(data_info_bit,Nsymbolperframe,bitPerSymbol));
        x = qammod(data_temp,QAM_mod,'gray');
        X = reshape(x,M,N);    % original QAM symbol mapped to Delay-Doppler Domain
        %% ISFFT
        X_TF = ISFFT(X);
        %% Heisenberg Transform
        X_ht = ifft(X_TF,M); 
        s = reshape(X_ht,[],1);
        
        %% Through Delay Doppler Channel
        r = zeros(size(s));
        for iTarget = 1:length(targetRange)
            temp = s.*(d.^lp(iTarget));
            r = r + (randn+1j*randn)*[temp(M*N-kp(iTarget)+1:M*N);temp(1:M*N-kp(iTarget))]/sqrt(2);
        end
        r = r + (randn(size(r)) + 1j*randn(size(r)))*sqrt(sigmaSquare(ii)/2);
        %% Wigner Transform
        Y = reshape(r,M,N);
        y_TF =  fft(Y,M);
        %% Sensing Based on TF domain 
        H_tf = y_TF .* conj(X_TF);
        rdm_tf = fft(ifft(H_tf).',N*10)'*sqrt(M/N);
        MM = max(abs(rdm_tf),[],'all');
        [I1,I2] = find(abs(rdm_tf)==MM);
        rangeEst = (I1-1)/(M*deltaf)*c/2;
        veloEst = (I2-1)/(N*10*T)*lambda/2;
        errorRange_TF(ii) = errorRange_TF(ii) + (rangeEst - targetRange)^2/NFrame;
        errorVelo_TF(ii)  = errorVelo_TF(ii)  + (veloEst  - targetSpeed)^2/NFrame;
% %         surf(abs(rdm_tf));
        %% SFFT
        y = SFFT(y_TF);
        y_vec = reshape(y,[],1);
        h_vec = ifft(conj(fft(x)).*fft(y_vec));
        H_est = reshape(h_vec,M,N);
        MM = max(abs(H_est),[],'all');
        [I1,I2] = find(abs(H_est)==MM);
        rangeEst = (I1-1)/(M*deltaf)*c/2;
        veloEst = (I2-1)/(N*T)*lambda/2;
        errorRange_MF(ii) = errorRange_MF(ii) + (rangeEst - targetRange)^2/NFrame;
        errorVelo_MF(ii)  = errorVelo_MF(ii)  + (veloEst  - targetSpeed)^2/NFrame;        
    end
end
errorRange_TF = sqrt(errorRange_TF);
errorVelo_TF = sqrt(errorVelo_TF);
errorRange_MF = sqrt(errorRange_MF);
errorVelo_MF = sqrt(errorVelo_MF);
figure(1);
semilogy(SNR_dB,errorRange_TF);
hold on;
semilogy(SNR_dB,errorRange_MF);
xlabel('SNR(dB)');
ylabel('Range RMSE(m)');
legend('Time Frequency Based Method','Match Filter Based Method');
savefig('fig/range_rmse_comparison.fig');
