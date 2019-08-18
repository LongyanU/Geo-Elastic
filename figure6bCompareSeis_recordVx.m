% first run figure5aStaggeredFDTra.m, figure5bStaggeredFD.m to get seismic records,
% then run figure6acompareSeimogramsVx,figure6bCompareSeis_recordVx; figure6cCompareSeis_recordTxx,
% figure6dcompareSeimogramsTxx to get the figures in figure 5 and figure 6. 
% This is only for the convenient of the reviewers.
clear;

clc
close all

load('figure5aStaggeredGridFDTra.mat')
figure;imagesc(-seis_recordVx(1:1600,45:end-45),[-0.002 0.002]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')



load('figure5bStaggeredGridFD.mat')
figure;imagesc(-seis_recordVx(1:1600,45:end-45),[-0.002 0.002]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


figure;imagesc(vs(45:end-45,45:end-45))
h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(xs-45,zs-45,'*r')
xlabel('x/dx')
ylabel('z/dz')
title('')

figure;imagesc(vp(45:end-45,45:end-45))
h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(xs-45,zs-45,'*r')
xlabel('x/dx')
ylabel('z/dz')
