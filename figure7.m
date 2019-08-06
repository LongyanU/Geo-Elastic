clear;
clc
close all
load('model2_VP_201X241_5m.mat')
figure;imagesc(cc)
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
set(get(h,'title'),'string','m/s');

load('model2_VS_201X241_5m.mat')
figure;imagesc(cc)
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
set(get(h,'title'),'string','m/s');