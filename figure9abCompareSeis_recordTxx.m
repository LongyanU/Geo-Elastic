% first run figure8aStaggeredFDTra.m, figure8bStaggeredFD.m to get seismic records,
% then run figure8cCompareSeis_recordVx, figure8dcompareSeimogramsVx; figure9abCompareSeis_recordTxx,
% figure9cdCompareSeimogramsTxx to get the figures in figure 8 and figure 9.
clear;
clc
close all


load('figure8aStaggeredGridFDTra.mat')
figure;imagesc(-seis_recordTxx(1:1600,45:end-45),[-2*10^3 2*10^3]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


load('figure8bStaggeredGridFD.mat')
figure;imagesc(-seis_recordTxx(1:1600,45:end-45),[-2*10^3 2*10^3]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

