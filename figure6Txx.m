
% First run figure5a.m, figure5b.m to get the snapshots and seismic records
% Then run Figure5Vx.m figure6Txx.m to plot the results used in the paper.
clear;
clc
close all

load TraLayered.mat
figure;imagesc(-seis_recordTxx(:,45:end-45),[-4*10^2 4*10^2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')
tt=seis_recordTxx;

% figure;plot(seis_recordTxx(:,65))
% xlabel('z/dz')
% ylabel('Amp')
% title('')
% axis([0 1600 -14*10^-3 17*10^-3])


load NewLayered.mat
figure;imagesc(-seis_recordTxx(:,45:end-45),[-4*10^2 4*10^2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


figure;plot(tt(:,65),'k-')
hold on; plot(seis_recordTxx(:,65),'r-')

title('')
xlabel('Travel time')
ylabel('Txx (Pa)')
legend('Tra FD scheme','Non-balanced FD scheme');
% axis([0 1600 -14*10^-3 24*10^-3])

% clip level 4