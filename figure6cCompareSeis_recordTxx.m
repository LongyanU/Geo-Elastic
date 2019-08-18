% first run figure5aStaggeredFDTra.m, figure5bStaggeredFD.m to get seismic records,
% then run figure6acompareSeimogramsVx,figure6bCompareSeis_recordVx; figure6cCompareSeis_recordTxx,
% figure6dcompareSeimogramsTxx to get the figures in figure 5 and figure 6. 
% This is only for the convenient of the reviewers.
clear;

clc
close all


load('figure5aStaggeredGridFDTra.mat')
figure;imagesc(-seis_recordTxx(1:1600,45:end-45),[-2*10^6 2*10^6]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


load('figure5bStaggeredGridFD.mat')
figure;imagesc(-seis_recordTxx(1:1600,45:end-45),[-2*10^6 2*10^6]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')



