clear;

clc
% close all

load TraLayered.mat
figure;imagesc(Txx);
% bb=Vx;
xlabel('x/dx')
ylabel('z/dz')
title('')

load NewLayered.mat
figure;imagesc(Txx);
xlabel('x/dx')
ylabel('z/dz')
title('')

figure;imagesc(vs)
h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(xs,zs,'*r')
xlabel('x/dx')
ylabel('z/dz')
title('')



% figure;imagesc(aa,[-10^-1.5 10^-1.5]);colormap gray;
% figure;imagesc(bb,[-10^-1.5 10^-1.5]);colormap gray;
% figure;imagesc(aa-bb,[-10^-1.5 10^-1.5]);colormap gray;
% 
% figure;plot(aa(:,70));
% hold on;plot(bb(:,70),'r');