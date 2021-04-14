clear all
close all
clc

M = 128;
data = 0:M-1;
sym = qammod(data,M,'gray');

scatterplot(sym,1,0,'b*');
for k = 1:M
    text(real(sym(k))-0.4,imag(sym(k))+0.4,dec2base(data(k),2,7));
end
axis([-12 12 -12 12])