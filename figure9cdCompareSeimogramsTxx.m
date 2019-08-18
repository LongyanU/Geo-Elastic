% first run figure8aStaggeredFDTra.m, figure8bStaggeredFD.m to get seismic records,
% then run figure8cCompareSeis_recordVx, figure8dcompareSeimogramsVx; figure9abCompareSeis_recordTxx,
% figure9cdCompareSeimogramsTxx to get the figures in figure 8 and figure 9.
clear;
clc
close all

load('figure8aStaggeredGridFDTra.mat')
figure;plot(seis_recordTxx(1:1600,65))
tt=seis_recordTxx;

load('figure8bStaggeredGridFD.mat')
hold  on;plot(seis_recordTxx(1:1600,65),'k')
tt2=seis_recordTxx;

tt3=tt2-tt;
hold on;plot(sign(real(tt3(1:1600,65))) .*sqrt(real(tt3(1:1600,65)).^2+ imag(tt3(1:1600,65)).^2) -0.3*10^4 ,'m')  %-0.2*10^-3

legend('Tra FD scheme','Non-balanced FD scheme', 'The difference between the  2 schemes');
title('')
grid on
xlabel('Travel time(ms)')
ylabel('Txx(Pa)')