
% First run figure8a.m, figure8b.m to get the snapshots and seismic records
% Then run Figure8Vx.m figure9Txx.m to plot the results used in the paper.
clear;

clc
close all

load('figure11a.mat')
figure;imagesc(-seis_recordTxx(:,45:end-45),[-4*10^-1 4*10^-1]);
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

load('figure11b.mat')
figure;imagesc(-seis_recordTxx(:,45:end-45),[-4*10^-1 4*10^-1]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

figure;plot(ttt(:,200),'k-')
hold on; plot(seis_recordTxx(:,200),'r-')
title('')
xlabel('Travel time')
ylabel('Txx (Pa)')
legend('Tra FD scheme','Non-balanced FD scheme');
axis([0 4500 -3 4.5])