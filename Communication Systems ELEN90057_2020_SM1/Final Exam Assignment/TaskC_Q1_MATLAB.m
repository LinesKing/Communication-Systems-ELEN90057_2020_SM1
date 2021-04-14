%% %%Task C: Multicarrier Modulation Q1%% %%

clear all
close all
clc

%% Input parameters
M = 16;
EsTx = 5; %% average energy per transimiting symbol EsTx = (M-1)/3*Eg
N0 = 0.2; %% power spectral density of the AWGN 
H = (1+1i)/sqrt(2); %% attenuation and phase shift

Eg = EsTx/(M-1)*3; %% Eg
d = sqrt(2*Eg); %% minimum distance
EsRx = (abs(H))^2*EsTx; %% average energy per received symbol

%% Random input
N = 100000; %% numbers of input
input_de = randi([0 M-1],N,1); %% decimal input
input_bi = de2bi(input_de); %% binary input

%% Gray-coded 16-QAM transmitter
QAM_bits = input_bi;
XR = zeros(1,size(QAM_bits,1)); 
XI = zeros(1,size(QAM_bits,1));
X = zeros(1,size(QAM_bits,1)); %% input sequence X

XR(and(QAM_bits(:,1) == 0,QAM_bits(:,2) == 0)) = -3/2*d;
XR(and(QAM_bits(:,1) == 0,QAM_bits(:,2) == 1)) = -1/2*d;
XR(and(QAM_bits(:,1) == 1,QAM_bits(:,2) == 1)) = 1/2*d;
XR(and(QAM_bits(:,1) == 1,QAM_bits(:,2) == 0)) = 3/2*d;

XI(and(QAM_bits(:,3) == 0,QAM_bits(:,4) == 0)) = 3/2*d;
XI(and(QAM_bits(:,3) == 0,QAM_bits(:,4) == 1)) = 1/2*d;
XI(and(QAM_bits(:,3) == 1,QAM_bits(:,4) == 1)) = -1/2*d;
XI(and(QAM_bits(:,3) == 1,QAM_bits(:,4) == 0)) = -3/2*d;

X = (XR + XI*1i); %% IQ mapping

%% IQ modulation
L = 16; %% symbol period
fc = 5; %% carrier frequency
fs = L * fc;    %% sampling frequency
t=0: 1/fs : (length(X)*fs-1)/fs; %% time vector

XL = upsample(X,fs); %% X(l) by upsampling X
F = ones(1,fs); %% pulse shaping
XL = filter(F,1,XL); %% transmitted signal X(l)
XLR = real(XL); %% real part of X(l)
XLI = imag(XL); %% imaginary part of X(l)

XLcR = XLR .* cos(2*pi*fc*t)/sqrt(fs/2); %% bandpass and normalization
XLcI = XLI .* (-sin(2*pi*fc*t))/sqrt(fs/2); %% bandpass and normalization

XLc = XLcR + XLcI*1i; %% transmitted signal with carrier

%% Channel transmission
H = H; %% attenuation and phase shift
mu = 0;
sigma = sqrt(N0/2); %% standard deviation of AWGN
WR = normrnd(mu,sigma,1,length(XLc)); %% real part of AWGN
WI = normrnd(mu,sigma,1,length(XLc)); %% imag part of AWGN
WL = WR + WI*1i;

YLc = XLc*H + WL; %% normalization of transmit power to one before AWGN

%% IQ demodulation
ZLc = YLc/H; %% Undo attenuation, phase shift and normalization
ZLcR = real(ZLc); %% real part of Z(l)
ZLcI = imag(ZLc); %% imaginary part of Z(l)

ZLcR_unmatched = sqrt(2)*ZLcR .* cos(2*pi*fc*t); %% sift and normalise
ZLcI_unmatched = sqrt(2)*ZLcI .* (-sin(2*pi*fc*t)); %% sift and normalise

F = ones(1,fs); %% matched filter
ZLR_matched = filter(fliplr(F),1,ZLcR_unmatched);
ZLI_matched = filter(fliplr(F),1,ZLcI_unmatched);

ZR = downsample(ZLR_matched,fs,fs-1)/sqrt(fs); %% downsample and normalise
ZI = downsample(ZLI_matched,fs,fs-1)/sqrt(fs); %% downsample and normalise

% constellation diagram and 16-QAM reference output (Gray-coded)
figure(1)
x = (0:15);
symgray = qammod(x,M,'gray');
scatter(ZR,ZI,5,'filled'); 

hold on
scatter(real(symgray)/sqrt(2),imag(symgray)/sqrt(2),'*');
% for k = 1:M
%     text(real(symgray(k)) - 0.0,imag(symgray(k)) + 0.3,...
%         dec2base(x(k),2,4),'Color',[1 0 0]);
%      text(real(symgray(k)) - 0.5,imag(symgray(k)) + 0.3,...
%          num2str(x(k)),'Color',[1 0 0]);
% end
grid on
xlabel('Inphase')
ylabel('Quadrature')
title('16-QAM constellation')

%% Decode 16-QAM transmitter
%Dicision
ZR(ZR < -d) = -3/2*d;
ZR(and(ZR > -d,ZR < 0)) = -1/2*d;
ZR(and(ZR > 0,ZR < d)) = 1/2*d;
ZR(ZR > d) = 3/2*d;

ZI(ZI < -d) = -3/2*d;
ZI(and(ZI > -d,ZI < 0)) = -1/2*d;
ZI(and(ZI > 0,ZI < d)) = 1/2*d;
ZI(ZI > d) = 3/2*d;

%Decoder
recv_bits = zeros(size(QAM_bits));

recv_bits(ZR < 0,1) = 0;
recv_bits(ZR > 0,1) = 1;
recv_bits(or(ZR == -3/2*d,ZR == 3/2*d),2) = 0;
recv_bits(or(ZR == -1/2*d,ZR == 1/2*d),2) = 1;
recv_bits(ZI > 0,3) = 0;
recv_bits(ZI < 0,3) = 1;
recv_bits(or(ZI == -3/2*d,ZI == 3/2*d),4) = 0;
recv_bits(or(ZI == -1/2*d,ZI == 1/2*d),4) = 1;

%% Output
output_bi = recv_bits;
output_de = bi2de(output_bi);
SER_simulated = 1 - mean(output_de == input_de) %% simulated SER

SER_sqrtM_PAM = 2*(1-1/sqrt(M))*qfunc(sqrt(3*EsRx/(M-1)/N0));
SER_theoretical = 1 - (1-SER_sqrtM_PAM)^2 %% theoretical SER
% SER_upbound = 4*(1-1/sqrt(M))*qfunc(sqrt(3*EsRx/(M-1)/N0))

