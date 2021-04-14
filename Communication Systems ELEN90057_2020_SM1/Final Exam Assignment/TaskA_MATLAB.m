%% %%Task A: 3-ary PAM%% %%

clear all
close all
clc

%% Q2
M = 3; 
a = 1;
b = 0.75;
Es = 2*a^2/3; %% average symbol energy
N = 100000; %% trials
data = randi([-(M-1)/2 (M-1)/2],N,1); %% s1=-1, s2=0, s3=1

EsN0dB = -4:20; %% Es/N0 in dB
EsN0 = 10.^(EsN0dB/20); %% linear Es/N0
N0 = Es./EsN0; 
SE_posterior = zeros(size(EsN0)); %% posterior symbol error
SEP_theoretical = zeros(size(EsN0)); %% theoretical symbol error probability

for i = 1:length(EsN0dB) 
    
    %calculate SE_posterior
    noise = normrnd(0,sqrt(N0(i)/2),size(data));
    for j = 1:length(data)
        switch data(j)
        case -1
            SE_posterior(i) = SE_posterior(i) + (noise(j) > (a-b));
        case 0
            SE_posterior(i) = SE_posterior(i) + ((noise(j) < (-b))||(noise(j) > (b)));
        case 1
            SE_posterior(i) = SE_posterior(i) + (noise(j) < (b-a));
        end
    end
    
    %calculate SEP_theoretical
    SEP_theoretical(i) = 2/3*qfunc(b/sqrt(N0(i)/2))+2/3*qfunc((a-b)/sqrt(N0(i)/2));
end
SEP_posterior = SE_posterior/N;%calculate SEP_posterior

figure
plot(EsN0dB,SEP_posterior,EsN0dB,SEP_theoretical);
grid on
legend('simulated SEP','theoretical SEP')
xlabel('Es/N0 / dB')
ylabel('SEP')
title('Symbol error probability curve for 3-PAM modulation')

%% Q4
M = 3; 
a = 1;
N0 = 0.1;
N = 100000;
b = 0:0.01:1;

%calculate SEP_theoretical
SEP_theoretical = 2/3*qfunc(b/sqrt(N0/2))+2/3*qfunc((a-b)/sqrt(N0/2));

figure
plot(b,SEP_theoretical)
grid on
xlabel('b')
ylabel('SEP')
title('Symbol error probability curve for 3-PAM modulation')

%% Q6
M = 3; 
a = 1;
N0 = 0.1;
N = 100000;
p = 0.01:0.01:0.49;
b = a/2-N0/2/a*log(p./(1-2.*p));

%calculate SEP_theoretical
SEP_theoretical = 2*(1-2*p).*qfunc(b/sqrt(N0/2))+2*p.*qfunc((a-b)/sqrt(N0/2));

figure
plot(p,SEP_theoretical)
grid on
xlabel('p')
ylabel('SEP')
title('Symbol error probability curve for 3-PAM modulation')





