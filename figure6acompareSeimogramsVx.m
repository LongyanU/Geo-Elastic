% first run figure5aStaggeredFDTra.m, figure5bStaggeredFD.m to get seismic records,
% then run figure6acompareSeimogramsVx,figure6bCompareSeis_recordVx; figure6cCompareSeis_recordTxx,
% figure6dcompareSeimogramsTxx to get the figures in figure 5 and figure 6. 
% This is only for the convenient of the reviewers.
clear;

clc
close all

load('figure5aStaggeredGridFDTra.mat')
figure;plot(seis_recordVx(1:1600,65))
tt=seis_recordVx;

load('figure5bStaggeredGridFD.mat')
hold  on;plot(seis_recordVx(1:1600,65),'k')
tt2=seis_recordVx;




tt3=tt2-tt;
hold on;plot(sign(real(tt3(1:1600,65))) .*sqrt(real(tt3(1:1600,65)).^2+ imag(tt3(1:1600,65)).^2) -0.001 ,'m')  %-0.2*10^-3

legend('Tra FD scheme','Non-balanced FD scheme', 'The difference between the 2 schemes');
title('')
grid on
xlabel('Travel time(ms)')
ylabel('Velocity (m/s)')
