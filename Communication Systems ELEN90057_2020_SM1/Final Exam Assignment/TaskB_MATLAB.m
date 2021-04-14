clear all
close all
clc

%% image encoder
img_data = imread('University-Welcome-sign.jpg');
dim_size = size(img_data);

img_bits = reshape(img_data,prod(dim_size),1);
img_bits = de2bi(img_bits);
img_bits = reshape(img_bits,1,prod(dim_size)*8);

%% Gray-coded 4-PAM transmitter
PAM_bits = double(img_bits);
PAM_bits = reshape(img_bits,[length(img_bits)/2,2]);

PAM_sym = zeros(1,size(PAM_bits,1));

PAM_sym(and(PAM_bits(:,1) == 0,PAM_bits(:,2) == 0)) = -3;
PAM_sym(and(PAM_bits(:,1) == 0,PAM_bits(:,2) == 1)) = -1;
PAM_sym(and(PAM_bits(:,1) == 1,PAM_bits(:,2) == 1)) = 1;
PAM_sym(and(PAM_bits(:,1) == 1,PAM_bits(:,2) == 0)) = 3;

%% channel transmission
channel_taps = [1,-.3,.45,-.25];

recv_signal = conv(PAM_sym,channel_taps);

%% %% %% Linear Equalizer %% %% %% 
equalized_signal = zeros(1,length(PAM_sym));

for L = 1:10
    h = zeros(L);
    hv = [channel_taps zeros(1,(L-4)*(L>4))]; %% h0 to hL
    
    for i = 0:L-1
        hdiag = zeros(L-i,1);
        hdiag(1:L-i,1) = hv(i+1); %% h diagonal
        h = h+diag(hdiag,-i);
    end
    c(1:L,L) = h\[1;zeros(L-1,1)]; %% c = h\q
end
display(c)

L = 4; %% L-tap ZF equalizer
equa_signal = conv(recv_signal,c(:,L)); %% equalizer
equalized_signal(L,1:length(PAM_sym)) = equa_signal(1:length(PAM_sym));
recv_signal = equalized_signal(L,1:length(PAM_sym));

%%%%%%%%%%%%%%%%%%%%%%
%% %% %% 4-PAM receiver (decision boundaries at -2,0,+2) %% %% %%
%Dicision
recv_signal(recv_signal<-2) = -3;
recv_signal((-2<recv_signal)&(recv_signal<0)) = -1;
recv_signal((0<recv_signal)&(recv_signal<2)) = 1;
recv_signal(recv_signal>2) = 3;

recv_bits = uint8(zeros(length(img_bits)/2,2));
%Decoder
recv_bits(or(recv_signal(:) == -3,recv_signal(:) == -1),1) = 0;
recv_bits(or(recv_signal(:) == 1,recv_signal(:) == 3),1) = 1;
recv_bits(or(recv_signal(:) == -3,recv_signal(:) == 3),2) = 0;
recv_bits(or(recv_signal(:) == -1,recv_signal(:) == 1),2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image decoder
recv_bits = reshape(recv_bits,size(img_bits));
recv_img = reshape(recv_bits,prod(dim_size),8);
recv_img = bi2de(recv_img);
recv_img = reshape(recv_img,dim_size);

figure(1);
imshow(recv_img)

%% Q2
SE = 0; %% symbol error

for Ik1 = -3:2:3
    for Ik2 = -3:2:3
        for Ik3 = -3:2:3
            if 0.3*Ik1-0.45*Ik2+0.25*Ik3 > 1
                SE = SE + 1;
            end
        end
    end
end

SEP = SE/64; %% symbol error probability

%% Q4
SE = 0; %% symbol error

for Ik1 = -3:2:3
    for Ik2 = -3:2:3
        for Ik3 = -3:2:3
            if 0.2391*Ik1-0.0932*Ik2+0.0018*Ik3 > 1
                SE = SE + 1;
            end
        end
    end
end

SEP = SE/64; %% symbol error probability

%% Q5
Rcha = recv_img(:,:,1) == img_data(:,:,1); %% Rchannel
Gcha = recv_img(:,:,2) == img_data(:,:,2); %% Gchannel
Bcha = recv_img(:,:,3) == img_data(:,:,3); %% Bchannel
ratio_err = 1 - mean((Rcha+Gcha+Bcha) == 3,'all');

figure(2)
plot(1:10,[0.9439 0.8206 0.4375 0.4288 zeros(1,6)])
grid on
xlabel('L')
ylabel('SEP')
title('Ratio of erroneous pixels and total number of pixels')

%% Q6
q = zeros(10,13);
for L = 1:10
    q(L,:) = conv(channel_taps,c(:,L)); %% q
end

SE = zeros(1,10); %% symbol error
for i=1:10
    for Ik1 = -3:2:3
        for Ik2 = -3:2:3
            for Ik3 = -3:2:3
                if  q(i,i+1)*Ik1+q(i,i+2)*Ik2+q(i,i+3)*Ik3 > 1
                    SE(i) = SE(i) + 1;
                end
            end
        end
    end
end

SEPtotal = 3/2*SE/64; %% symbol error probability

figure(3)
plot(1:10,SEPtotal)
grid on
xlabel('L')
ylabel('SEPtotal')
title('Exact symbol error probability against the equaliser length L')
