%% %%Task C: Multicarrier Modulation Q2%% %%
clear all
close all
clc

%% Input parameters
N = 10; %% N independent sub-channels
n = 1:N;
Pe = 1e-2; %% the highest acceptable symbol-error-rate
EsTxN0dB = 60; %% EsTx/N0 in dB
H2dB = -30-3*n; %% square of attenuation on each sub-channel

EsTxN0linear = 10.^(EsTxN0dB/10); %% linear EsTx/N0
Hlinear = sqrt(10.^(H2dB/10)); %% linear attenuation on each sub-channel
EsRxN0linear = (abs(Hlinear)).^2*EsTxN0linear; %% average energy per received symbol

% find Mmax
M = 64;
SER_sqrtM_PAM = 2*(1-1/sqrt(M))*qfunc(sqrt(3.*EsRxN0linear/(M-1)));
SER_theoretical = 1 - (1-SER_sqrtM_PAM).^2

Mmax = 64; %% calculate maxima of Mn
bmax = log2(Mmax); %% calculate maxima of bn
brange = bmax:-2:0; %% range of bn
Mrange = 2.^(brange); %% range of Mn

% find Mmax
for i = 1:N
    M1(i) = QAM_sub_channel(EsTxN0linear,Hlinear(i),Pe,Mrange);
end

channel_use = M1>1; %% find used channel
display(channel_use)
EsTxN0reallocate = EsTxN0linear*N/sum(channel_use); %% reallocate the energy equally between the remaining sub-channels
EsRxN0reallocate = (abs(Hlinear)).^2*EsTxN0reallocate; %% average energy per received symbol

% find Mmax for reallocation
M = 256;
SER_sqrtM_PAM = 2*(1-1/sqrt(M))*qfunc(sqrt(3.*EsRxN0reallocate/(M-1)));
SER_theoretical = 1 - (1-SER_sqrtM_PAM).^2

% find Mmax again
M2(1) = M;
for i = 2:sum(channel_use)
    M2(i) = QAM_sub_channel(EsTxN0reallocate,Hlinear(i),Pe,Mrange);
end

% find the bits tansmitted 
bits_tansmitted1 = sum(log2(M1));
bits_tansmitted2 = sum(log2(M2));
display([M1 bits_tansmitted1])
display([M2 bits_tansmitted2])