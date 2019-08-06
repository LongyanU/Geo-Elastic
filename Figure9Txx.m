
% First run figure8a.m, figure8b.m to get the snapshots and seismic records
% Then run Figure8Vx.m figure9Txx.m to plot the results used in the paper.
clear;

clc
close all

load('figure8a.mat')
figure;imagesc(-seis_recordTxx(:,45:end-45),[-4*10^2 4*10^2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

% figure;plot(seis_recordTxx(:,65))
% xlabel('z/dz')
% ylabel('Amp')
% title('')
% axis([0 1600 -14*10^-3 17*10^-3])
ttt=seis_recordTxx;

load('figure8b.mat')
figure;imagesc(-seis_recordTxx(:,45:end-45),[-4*10^2 4*10^2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

figure;plot(ttt(:,65),'k-')
hold on; plot(seis_recordTxx(:,65),'r-')
title('')
xlabel('Travel time')
ylabel('Txx (Pa)')
legend('Tra FD scheme','Non-balanced FD scheme');
axis([0 2000 -1*10^4 1.6*10^4])