%% %%Task C: Multicarrier Modulation Q3%% %%
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

% % find Mmax
% M = 128;
% SER_sqrtM_PAM = 2*(1-1/sqrt(M))*qfunc(sqrt(3.*EsRxN0linear/(M-1)));
% SER_theoretical = 1 - (1-SER_sqrtM_PAM).^2

Mmax = 128; %% calculate maxima of Mn
bmax = log2(Mmax); %% calculate maxima of bn
brange = bmax:-1:0; %% range of bn
Mrange = 2.^(brange); %% range of Mn

% find Mmax
for i = 1:N
    M1(i) = QAM_sub_channel2(EsTxN0linear,Hlinear(i),Pe,Mrange);
end

channel_use = M1>1; %% find used channel
display(channel_use)
EsRxN0reallocate = EsTxN0linear*N/sum(channel_use); %% reallocate the energy equally between the remaining sub-channels
EsRxN0reallocatedB = 10*log10(EsRxN0reallocate);
% find Mmax again
for i = 1:sum(channel_use)
    M2(i) = QAM_sub_channel2(EsRxN0reallocate,Hlinear(i),Pe,Mrange);
end

% % find the minimum energy E for sub-channel 1 to operate 128-ary-QAM
% for i = 1:bmax
%     for j = 1:N
%         c = (1-sqrt(1-Pe))/2/(1-1/sqrt(Mrange(i)));
%         EsTxN0linear1 = qfuncinv(c)^2*(Mrange(i)-1)/3/(Hlinear(j))^2;
%         Emin(j,i) = 10*log10(EsTxN0linear1);
%     end
% end

% new reallocation strategy
EsRxN0reallocate2dB = [58.2 61.2 61.3 61 60.9 61.6 59.6 59.8];
EsRxN0reallocate2 = 10.^(EsRxN0reallocate2dB/10);
flag = sum(EsRxN0reallocate2) < 10e6
for i = 1:length(EsRxN0reallocate2dB)
    M3(i) = QAM_sub_channel2(EsRxN0reallocate2(i),Hlinear(i),Pe,Mrange);
end

% find the bits tansmitted 
bits_tansmitted1 = sum(log2(M1));
bits_tansmitted2 = sum(log2(M2));
bits_tansmitted3 = sum(log2(M3));
display([M1 bits_tansmitted1])
display([M2 bits_tansmitted2])
display([M3 bits_tansmitted3])


