% first run Figure11aStaggeredFDTra.m,figure11bStaggeredFD.m,figure11cKspace to get seismic records,
% then run figure11compareSeimogramsVz.m, Figure11CompareSeis_recordVz.m; figure12compareSeimogramsTxx.m,figure12CompareSeis_recordTxx
% to get the figures in figure 11 and figure 12.
% This is only for the convenient of the reviewers.

clear;

clc
close all


load('StaggeredGridFDTra2.mat')
figure;imagesc(-seis_recordTxx(1:4000,45:end-45),[-0.2 0.2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')



load('StaggeredGridFD2.mat')
figure;imagesc(-seis_recordTxx(1:4000,45:end-45),[-0.2 0.2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

load('KSpaceLayer2.mat')
figure;imagesc(-sign(real(seis_recordTxx(1:4000,45:end-45))).*sqrt(real(seis_recordTxx(1:4000,45:end-45)).^2+ imag(seis_recordTxx(1:4000,45:end-45)).^2 ), [-0.2 0.2]);
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')

