% first run Figure11aStaggeredFDTra.m,figure11bStaggeredFD.m,figure11cKspace to get seismic records,
% then run figure11compareSeimogramsVz.m, Figure11CompareSeis_recordVz.m; figure12compareSeimogramsTxx.m,figure12CompareSeis_recordTxx
% to get the figures in figure 11 and figure 12.
% This is only for the convenient of the reviewers.

clear;

clc
close all

load('StaggeredGridFDTra2.mat')
figure;plot(seis_recordTxx(1:4000,400))
tt=seis_recordTxx;

load('StaggeredGridFD2.mat')
hold  on;plot(seis_recordTxx(1:4000,400),'k')
tt2=seis_recordTxx;


load('KSpaceLayer2.mat')
hold on;plot(sign(real(seis_recordTxx(1:4000,400))) .*sqrt(real(seis_recordTxx(1:4000,400)).^2+ imag(seis_recordTxx(1:4000,400)).^2)   ,'r')


tt3=tt2-tt;
hold on;plot(sign(real(tt3(1:4000,400))) .*sqrt(real(tt3(1:4000,400)).^2+ imag(tt3(1:4000,400)).^2) -0.2 ,'m')  %-0.2*10^-3

legend('Tra FD scheme','Non-balanced FD scheme','K space', 'The difference between the first 2 schemes');
title('')
grid on
% xlabel('travel time');
xlabel('Travel time(ms)')
ylabel('Txx(Pa)')

figure;plot(sign(real(tt3(1:4000,400))) .*sqrt(real(tt3(1:4000,400)).^2+ imag(tt3(1:4000,400)).^2) ,'m')  %-0.2*10^-3