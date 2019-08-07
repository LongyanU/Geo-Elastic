
% First run figure8a.m, figure8b.m to get the snapshots and seismic records
% Then run Figure8Vx.m figure9Txx.m to plot the results used in the paper.
clear;

clc
close all

 load('figure11a.mat')
figure;imagesc(-seis_recordVx(:,45:end-45),[-0.0002 0.0002]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

% figure;plot(seis_recordVx(:,65))
% xlabel('z/dz')
% ylabel('Amp')
% title('')
% axis([0 2000 -7*10^-3 10.5*10^-3])

tt=seis_recordVx;


load('figure11b.mat')
figure;imagesc(-seis_recordVx(:,45:end-45),[-0.0002 0.0002]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


figure;plot(tt(:,200),'k-')
hold on; plot(seis_recordVx(:,200),'r-')
xlabel('Travel time(s)')
ylabel('Velocity (m/s)')
legend('Tra FD scheme','Non-balanced FD scheme');
title('')
axis([0 4500 -2*10^-3 3*10^-3])