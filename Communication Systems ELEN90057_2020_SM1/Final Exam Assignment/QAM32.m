function PER = QAM32(EsTx,N0,H,M)

    % Input parameters
    d = sqrt(EsTx/5); %% dmin
    
    % Random input
    N = 10000; %% numbers of input
    input_de = randi([0 M-1],N,1);
    input_bi = de2bi(input_de);
    
    % Gray-coded 32-QAM transmitter
    QAM_bits = input_bi;
    XR = zeros(1,size(QAM_bits,1));
    XI = zeros(1,size(QAM_bits,1));
    X = zeros(1,size(QAM_bits,1));

    XR((QAM_bits(:,1) == 0)&(QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 1)) = -5/2*d;
    XR((QAM_bits(:,1) == 0)&(QAM_bits(:,2) == 1)&(QAM_bits(:,3) == 1)) = -3/2*d;
    XR((QAM_bits(:,1) == 0)&(QAM_bits(:,2) == 1)&(QAM_bits(:,3) == 0)) = -1/2*d;
    XR((QAM_bits(:,1) == 1)&(QAM_bits(:,2) == 1)&(QAM_bits(:,3) == 0)) = 1/2*d;
    XR((QAM_bits(:,1) == 1)&(QAM_bits(:,2) == 1)&(QAM_bits(:,3) == 1)) = 3/2*d;
    XR((QAM_bits(:,1) == 1)&(QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 1)) = 5/2*d;
    
    XI((QAM_bits(:,4) == 0)&(QAM_bits(:,5) == 0)) = 3/2*d;
    XI((QAM_bits(:,4) == 0)&(QAM_bits(:,5) == 1)) = 1/2*d;
    XI((QAM_bits(:,4) == 1)&(QAM_bits(:,5) == 1)) = -1/2*d;
    XI((QAM_bits(:,4) == 1)&(QAM_bits(:,5) == 0)) = -3/2*d;
    
    XR((QAM_bits(:,1) == 0)&(QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 0)&(QAM_bits(:,5) == 0)) = -3/2*d;
    XR((QAM_bits(:,1) == 0)&(QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 0)&(QAM_bits(:,5) == 1)) = -1/2*d;
    XR((QAM_bits(:,1) == 1)&(QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 0)&(QAM_bits(:,5) == 1)) = 1/2*d;
    XR((QAM_bits(:,1) == 1)&(QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 0)&(QAM_bits(:,5) == 0)) = 3/2*d;
    
    XI((QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 0)&(QAM_bits(:,4) == 0)) = 5/2*d;
    XI((QAM_bits(:,2) == 0)&(QAM_bits(:,3) == 0)&(QAM_bits(:,4) == 1)) = -5/2*d;

    X = (XR + XI*1i); %% mapping
    
    % IQ modulation
    L = 16; %% symbol period
    fc = 5; %% carrier frequency
    fs = L * fc;    % sampling frequency
    t=0: 1/fs : (length(X)*fs-1)/fs; %time vector

    XL = upsample(X,fs); 
    F = ones(1,fs);
    XL = filter(F,1,XL); %% transmitted signal
    XLR = real(XL);
    XLI = imag(XL);

    XLR = XLR .* cos(2*pi*fc*t)/sqrt(fs/2);
    XLI = XLI .* (-sin(2*pi*fc*t))/sqrt(fs/2);

    XLc = XLR + XLI*1i; %% transmitted signal with carrier
    
    % Channel transmission
    mu = 0;
    sigma = sqrt(N0/2);
    WR = normrnd(mu,sigma,1,length(XLc)); %% real part noise
    WI = normrnd(mu,sigma,1,length(XLc)); %% imag part noise
    WL = WR + WI*1i;

    YLc = XLc*H + WL; %% normalization of transmit power to one
    
    % IQ demodulation
    ZLc = YLc/H; %% Undo attenuation and phase shift plus normalization
    ZLcR = real(ZLc);
    ZLcI = imag(ZLc);

    ZLcR_unmatched = sqrt(2)*ZLcR .* cos(2*pi*fc*t);
    ZLcI_unmatched = sqrt(2)*ZLcI .* (-sin(2*pi*fc*t));

    F = ones(1,fs); %% matched filter
    ZLR_matched = filter(fliplr(F),1,ZLcR_unmatched);
    ZLI_matched = filter(fliplr(F),1,ZLcI_unmatched);

    ZR = downsample(ZLR_matched,fs,fs-1)/sqrt(fs); %% sample real and normalise
    ZI = downsample(ZLI_matched,fs,fs-1)/sqrt(fs); %% sample imag and normalise

    % Decode 32-QAM transmitter
    ZR(ZR < -2*d) = -5/2*d;
    ZR(and(ZR > -2*d,ZR < -d)) = -3/2*d;
    ZR(and(ZR > -d,ZR < 0)) = -1/2*d;
    ZR(and(ZR > 0,ZR < d)) = 1/2*d;
    ZR(and(ZR > d,ZR < 2*d)) = 3/2*d;
    ZR(ZR > 2*d) = 5/2*d;

    ZI(ZI < -2*d) = -5/2*d;
    ZI(and(ZI > -2*d,ZI < -d)) = -3/2*d;
    ZI(and(ZI > -d,ZI < 0)) = -1/2*d;
    ZI(and(ZI > 0,ZI < d)) = 1/2*d;
    ZI(and(ZI > d,ZI < 2*d)) = 3/2*d;
    ZI(ZI > 2*d) = 5/2*d;

    recv_bits = zeros(size(QAM_bits));

    recv_bits(ZR < 0,1) = 0;
    recv_bits(ZR > 0,1) = 1;
    recv_bits(or(abs(ZR) == 5/2*d,abs(ZI) == 5/2*d),2) = 0;
    recv_bits(and(abs(ZR) ~= 5/2*d,abs(ZI) ~= 5/2*d),2) = 1;
    recv_bits(or(abs(ZR) == 1/2*d,abs(ZI) == 5/2*d),3) = 0;
    recv_bits(and(abs(ZR) ~= 1/2*d,abs(ZI) ~= 5/2*d),3) = 1;
    recv_bits(ZI > 0,4) = 0;
    recv_bits(ZI < 0,4) = 1;
    recv_bits((((abs(ZI) == 5/2*d)&(abs(ZR) == 3/2*d))|(abs(ZI) == 3/2*d)),5) = 0;
    recv_bits((((abs(ZI) == 5/2*d)&(abs(ZR) == 1/2*d))|(abs(ZI) == 1/2*d)),5) = 1;

    % Output
    output_bi = recv_bits;
    output_de = bi2de(output_bi);

    PER = 1 - mean(output_de == input_de);
end

