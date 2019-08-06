
% First run figure5a.m, figure5b.m to get the snapshots and seismic records
% Then run Figure5Vx.m figure6Txx.m to plot the results used in the paper.

clear;

clc
close all

load TraLayered.mat
figure;imagesc(-seis_recordVx(:,45:end-45),[-0.0002 0.0002]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

% figure;plot(seis_recordVx(:,65))
% xlabel('z/dz')
% ylabel('Amp')
% title('')
% % axis([0 1600 -12*10^-3 17*10^-3])
tt=seis_recordVx;

load NewLayered.mat
figure;imagesc(-seis_recordVx(:,45:end-45),[-0.0002 0.0002]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


figure;plot(tt(:,65),'k-')
hold on; plot(seis_recordVx(:,65),'r-')
xlabel('Travel time')
ylabel('Velocity (m/s)')
legend('Tra FD scheme','Non-balanced FD scheme');
title('')
% axis([0 1600 -14*10^-3 22*10^-3])