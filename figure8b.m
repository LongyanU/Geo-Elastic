% First run figure8a.m, figure8b.m to get the snapshots and seismic records
% Then run Figure8Vx.m figure9Txx.m to plot the results used in the paper.

clear;
clc
close all;

load('model2_VP_201X241_5m.mat')
vp=cc;
[nz nx]=size(cc);
nx=nx+45*2;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=vp(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=vp(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=vp(ii,201);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end
clear vp;
vp=vv;

load('model2_VS_201X241_5m.mat')
vs=cc;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=vs(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=vs(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=vs(ii,201);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end
clear vs;
vs=vv;



load('model2_DEN_201X241_5m.mat')
rou=cc;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=rou(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=rou(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=rou(ii,201);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end
clear rou;
rou=vv;

rou=rou*1000;



nt=2000;
isnap=10;    % snapshot sampling
h=10;
dx=h;


x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;

f0=78;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


% Source location
xs=floor(nz/2);
zs=46;

h=dx;

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics

tic

p=zeros([nz nx]); Vx=p; Vz=p;
Txxx=p;
Txzz=p;
Tzzz=p;
Txzx=p;

Txx=p;
Txz=p;
Tzz=p;


Vxx=p;
Vzz=p;
Vxz=p;
Vzx=p;

coeff=[1.57665 -0.290813 0.0879337 -0.0293731 0.00924336 -0.00239015 0.000409115];
% coeff=[ 1.25181, -0.121753, 0.0335093, -0.0110582, 0.00346284, -0.000854221, 0.000118311];
for it=1:nt-2,
    
    %Txx/x
    Txxx=coeff(1)*(Txx-circshift(Txx,[0 1]))+...
        coeff(2)*(circshift(Txx,[0 -1])-circshift(Txx,[0 2]))+...
        coeff(3)*(circshift(Txx,[0 -2])-circshift(Txx,[0 3]))+...
        coeff(4)*(circshift(Txx,[0 -3])-circshift(Txx,[0 4]))+...
        coeff(5)*(circshift(Txx,[0 -4])-circshift(Txx,[0 5]))+...
        coeff(6)*(circshift(Txx,[0 -5])-circshift(Txx,[0 6]))+...
        coeff(7)*(circshift(Txx,[0 -6])-circshift(Txx,[0 7]));
    
    %Txz/z
    Txzz=coeff(1)*(circshift(Txz,[-1 0])-Txz)+...
        coeff(2)*(circshift(Txz,[-2 0])-circshift(Txz,[1 0]))+...
        coeff(3)*(circshift(Txz,[-3 0])-circshift(Txz,[2 0]))+...
        coeff(4)*(circshift(Txz,[-4 0])-circshift(Txz,[3 0]))+...
        coeff(5)*(circshift(Txz,[-5 0])-circshift(Txz,[4 0]))+...
        coeff(6)*(circshift(Txz,[-6 0])-circshift(Txz,[5 0]))+...
        coeff(7)*(circshift(Txz,[-7 0])-circshift(Txz,[6 0]));
    
    %Tzz/z
    Tzzz=coeff(1)*(Tzz-circshift(Tzz,[1 0]))+...
        coeff(2)*(circshift(Tzz,[-1 0])-circshift(Tzz,[2 0]))+...
        coeff(3)*(circshift(Tzz,[-2 0])-circshift(Tzz,[3 0]))+...
        coeff(4)*(circshift(Tzz,[-3 0])-circshift(Tzz,[4 0]))+...
        coeff(5)*(circshift(Tzz,[-4 0])-circshift(Tzz,[5 0]))+...
        coeff(6)*(circshift(Tzz,[-5 0])-circshift(Tzz,[6 0]))+...
        coeff(7)*(circshift(Tzz,[-6 0])-circshift(Tzz,[7 0]));
    
    
    %Txz/x
    Txzx=coeff(1)*(circshift(Txz,[0 -1])-Txz)+...
        coeff(2)*(circshift(Txz,[0 -2])-circshift(Txz,[0 1]))+...
        coeff(3)*(circshift(Txz,[0 -3])-circshift(Txz,[0 2]))+...
        coeff(4)*(circshift(Txz,[0 -4])-circshift(Txz,[0 3]))+...
        coeff(5)*(circshift(Txz,[0 -5])-circshift(Txz,[0 4]))+...
        coeff(6)*(circshift(Txz,[0 -6])-circshift(Txz,[0 5]))+...
        coeff(7)*(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));
    
    
    Vx=Vx+1./(rou).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rou).*dt.*(Tzzz+Txzx)/h;
    Vx(zs,xs)=Vx(zs,xs)+src(it);
    Vz(zs,xs)=Vz(zs,xs)+src(it);
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
    Vxx=(circshift(Vx,[0 -1])-circshift(Vx,[0 0]));
    %         coeff(2)*(circshift(Vx,[0 -2])-circshift(Vx,[0 1]))+...
    %         coeff(3)*(circshift(Vx,[0 -3])-circshift(Vx,[0 2]))+...
    %         coeff(4)*(circshift(Vx,[0 -4])-circshift(Vx,[0 3]))+...
    %         coeff(5)*(circshift(Vx,[0 -5])-circshift(Vx,[0 4]))+...
    %         coeff(6)*(circshift(Vx,[0 -6])-circshift(Vx,[0 5]))+...
    %         coeff(7)*(circshift(Vx,[0 -7])-circshift(Vx,[0 6]));
    
    Vxz=(circshift(Vx,[0])-circshift(Vx,[1]));
    %         coeff(2)*(circshift(Vx,[-1])-circshift(Vx,[2]))+...
    %         coeff(3)*(circshift(Vx,[-2])-circshift(Vx,[3]))+...
    %         coeff(4)*(circshift(Vx,[-3])-circshift(Vx,[4]))+...
    %         coeff(5)*(circshift(Vx,[-4])-circshift(Vx,[5]))+...
    %         coeff(6)*(circshift(Vx,[-5])-circshift(Vx,[6]))+...
    %         coeff(7)*(circshift(Vx,[-6])-circshift(Vx,[7]));
    
    Vzz=(circshift(Vz,[-1])-circshift(Vz,[0]));
    %         coeff(2)*(circshift(Vz,[-2])-circshift(Vz,[1]))+...
    %         coeff(3)*(circshift(Vz,[-3])-circshift(Vz,[2]))+...
    %         coeff(4)*(circshift(Vz,[-4])-circshift(Vz,[3]))+...
    %         coeff(5)*(circshift(Vz,[-5])-circshift(Vz,[4]))+...
    %         coeff(6)*(circshift(Vz,[-6])-circshift(Vz,[5]))+...
    %         coeff(7)*(circshift(Vz,[-7])-circshift(Vz,[6]));
    
    
    Vzx=(circshift(Vz,[0 0])-circshift(Vz,[0 1]));
    %         coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
    %         coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
    %         coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
    %         coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
    %         coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
    %         coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    
    
%     Txx=Txx+dt*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
%     Tzz=Tzz+dt*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
%     Txz=Txz+dt*(vs.^2).*(Vxz+Vzx)/h;
    
     Txx=Txx+(rou).*dt.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+(rou).*dt.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+(rou).*dt.*(vs.^2).*(Vxz+Vzx)/h;
    
    seis_recordTxx(it,:)=Txx(zs,:);
    seis_recordTzz(it,:)=Tzz(zs,:);
    seis_recordTxz(it,:)=Txz(zs,:);
    
    
    if rem(it,isnap)== 0,
        imagesc(x,z,Vx), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
    
    
    
    
    %     if it==700
    %         pp1=Vx;
    %     elseif it==850
    %         pp2=Vx;
    %     end
    
     if it==1000
        pp1=Txx;
         Vx1=Vx;
    elseif it==1200
        pp2=Txx;
        Vx2=Vx;
    end
    
end


toc
%
save('figure8b.mat')