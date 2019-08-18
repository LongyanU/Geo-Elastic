%Equation 38's correctness
%When Vp== 5003, it is stable.
%When Vp=5005, it is unstable.
clear;
clc
close all;
nt=1602;
isnap=10;    % snapshot sampling
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
rou=ones(nz,nx);
f0=70;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;              % time derivative to obtain gaussian
 
for i=1:200
    for j=1:200
%         vp(i,j)=5003;  %stable
        vp(i,j)=5005;  %unstable
    end
end
 

xs=nz/2;
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
 
vs=vp./sqrt(3);
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
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
    Vxx=(circshift(Vx,[0 -1])-circshift(Vx,[0 0]));
    
    Vxz=(circshift(Vx,[0])-circshift(Vx,[1]));

    Vzz=(circshift(Vz,[-1])-circshift(Vz,[0]));
   
    Vzx=(circshift(Vz,[0 0])-circshift(Vz,[0 1]));
   
    
    Txx=Txx+dt*rou.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rou.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rou.*(vs.^2).*(Vxz+Vzx)/h;
    
 
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
