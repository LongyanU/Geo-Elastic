% first run figure5aStaggeredFDTra.m, figure5bStaggeredFD.m to get seismic records,
% then run figure6acompareSeimogramsVx,figure6bCompareSeis_recordVx; figure6cCompareSeis_recordTxx,
% figure6dcompareSeimogramsTxx to get the figures in figure 5 and figure 6. 
% This is only for the convenient of the reviewers.
% 时间已过 6.282245 秒。
clear;
clc
close all;
nt=1602;
isnap=100;    % snapshot sampling
h=10;
dx=h;
nx=200;
nz=200;

x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;

vp=zeros(nz,nx);

vs=vp;
rou=vp;
f0=78;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian



for i=1:100
    for j=1:200
        vp(i,j)=1500;
    end
end

for i=101:130
    for j=1:200
        vp(i,j)=2190;
    end
end

for i=131:135
    for j=1:200
        vp(i,j)=3126;
    end
end

for i=136:200
    for j=1:200
        vp(i,j)=2190;
    end
end

for i=101:130
    for j=1:200
        vs(i,j)=964;
    end
end

for i=131:135
    for j=1:200
        vs(i,j)=1851;
    end
end

for i=136:200
    for j=1:200
        vs(i,j)=964;
    end
end


for i=1:100
    for j=1:200
        rou(i,j)=1032;
    end
end

for i=101:130
    for j=1:200
        rou(i,j)=2085;
    end
end

for i=131:135
    for j=1:200
        rou(i,j)=2066;
    end
end

for i=136:200
    for j=1:200
        rou(i,j)=2085;
    end
end
% Source location
rou=rou*1000;
xs=nz/2;
zs=46;

h=dx;


seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordTxx=zeros(nt,nx);
seis_recordTzz=zeros(nt,nx);
seis_recordTxz=zeros(nt,nx);

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics

dx=h;
dz=h;


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


coeff=[ 1.25219, -0.121591, 0.0331761, -0.0109247, 0.0035117, -0.000949881, 0.000177978];
tic
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
    Txzz=coeff(1)*(circshift(Txz,[ 1])-circshift(Txz,[ 0]))+...
        coeff(2)*(circshift(Txz,[ 2])-circshift(Txz,[ -1]))+...
        coeff(3)*(circshift(Txz,[ 3])-circshift(Txz,[ -2]))+...
        coeff(4)*(circshift(Txz,[ 4])-circshift(Txz,[ -3]))+...
        coeff(5)*(circshift(Txz,[ 5])-circshift(Txz,[ -4]))+...
        coeff(6)*(circshift(Txz,[ 6])-circshift(Txz,[ -5]))+...
        coeff(7)*(circshift(Txz,[ 7])-circshift(Txz,[ -6]));
    
    %Tzz/z
    Tzzz=coeff(1)*(circshift(Tzz,[ 0])-circshift(Tzz,[ -1]))+...
        coeff(2)*(circshift(Tzz,[ 1])-circshift(Tzz,[ -2]))+...
        coeff(3)*(circshift(Tzz,[ 2])-circshift(Tzz,[ -3]))+...
        coeff(4)*(circshift(Tzz,[ 3])-circshift(Tzz,[ -4]))+...
        coeff(5)*(circshift(Tzz,[ 4])-circshift(Tzz,[ -5]))+...
        coeff(6)*(circshift(Tzz,[ 5])-circshift(Tzz,[ -6]))+...
        coeff(7)*(circshift(Tzz,[ 6])-circshift(Tzz,[ -7]));
    
    
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
    
   Vxx=coeff(1)*(circshift(Vx,[0 -1])-Vx)+...
        coeff(2)*(circshift(Vx,[0 -2])-circshift(Vx,[0 1]))+...
        coeff(3)*(circshift(Vx,[0 -3])-circshift(Vx,[0 2]))+...
        coeff(4)*(circshift(Vx,[0 -4])-circshift(Vx,[0 3]))+...
        coeff(5)*(circshift(Vx,[0 -5])-circshift(Vx,[0 4]))+...
        coeff(6)*(circshift(Vx,[0 -6])-circshift(Vx,[0 5]))+...
        coeff(7)*(circshift(Vx,[0 -7])-circshift(Vx,[0 6]));
    
    Vxz=coeff(1)*(circshift(Vx,[ 0])-circshift(Vx,[ -1]))+...
        coeff(2)*(circshift(Vx,[ 1])-circshift(Vx,[ -2]))+...
        coeff(3)*(circshift(Vx,[ 2])-circshift(Vx,[ -3]))+...
        coeff(4)*(circshift(Vx,[ 3])-circshift(Vx,[ -4]))+...
        coeff(5)*(circshift(Vx,[ 4])-circshift(Vx,[ -5]))+...
        coeff(6)*(circshift(Vx,[ 5])-circshift(Vx,[ -6]))+...
        coeff(7)*(circshift(Vx,[ 6])-circshift(Vx,[ -7]));
    
      Vzz=coeff(1)*(circshift(Vz,[ 1])-circshift(Vz,[ 0]))+...
        coeff(2)*(circshift(Vz,[ 2])-circshift(Vz,[ -1]))+...
        coeff(3)*(circshift(Vz,[ 3])-circshift(Vz,[ -2]))+...
        coeff(4)*(circshift(Vz,[ 4])-circshift(Vz,[ -3]))+...
        coeff(5)*(circshift(Vz,[ 5])-circshift(Vz,[ -4]))+...
        coeff(6)*(circshift(Vz,[ 6])-circshift(Vz,[ -5]))+...
        coeff(7)*(circshift(Vz,[ 7])-circshift(Vz,[ -6]));
    
    Vzx=coeff(1)*(Vz-circshift(Vz,[0 1]))+...
        coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
        coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
        coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
        coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
        coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
        coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    
    
    Txx=Txx+dt*rou.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rou.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rou.*(vs.^2).*(Vxz+Vzx)/h;
    
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
end

toc
save('figure5aStaggeredGridFDTra.mat')