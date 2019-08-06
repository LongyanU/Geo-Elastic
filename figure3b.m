clear;
clc

M=7;

AA=zeros(M,M);
b=zeros(M,1);
dt=0.001;

v=1500;
h=20;
k=linspace(1/50,0.83*pi/h,M);

for ii=1:M
    for kk=1:5
        xita=(kk-1)*pi/16;
        for jj=1:M
            AA(ii,jj)=2*cos(jj*k(ii)*h*cos(xita))-2*cos((jj-1)*k(ii)*h*cos(xita))...
                +2*cos(jj*k(ii)*h*sin(xita))-2*cos((jj-1)*k(ii)*h*sin(xita)) +AA(ii,jj);
        end
        b(ii)=1/(v^2*dt^2/h^2)*(2*cos(v*k(ii)*dt)-2)+b(ii);
    end
end
v=1500; %%figure 3a
digits(6)
c=AA\b
vpa(c)'
tau=dt;

v=2598;
r=v*tau/h;

k=linspace((pi)/(100*h),(pi)/h,100);
for kk=1:5
    xita=(kk-1)*pi/16;
    a=0;
    for n=1:M
        a=a+c(n)*sin((n-0.5)*k*h*cos(xita));
    end
    
    b=0;
    for n=1:M
        b=b+c(n)*sin((n-0.5)*k*h*sin(xita));
    end
    
    a=a.*sin((1-0.5)*k*h*cos(xita));
    b=b.*sin((1-0.5)*k*h*sin(xita));
    
    delta=2./(r*k*h);
    delta=delta.*asin(r*sqrt(a+b));
%     a1=(h/v*(1./delta-1));
    
    a1=delta;
    if kk==1  
        figure;plot(a1,'k','LineWidth',2)
    elseif kk==2
        hold on;plot(a1,'m--','LineWidth',2);
    elseif kk==3
        hold on;plot(a1,'r:','LineWidth',2)
    elseif kk==4
        hold on; plot(a1,'b-.','LineWidth',2)
    else
        hold on;plot(a1,'c:.','LineWidth',2)
    end
end
grid on
axis([0 100 0.97  1.01])
legend('¦È=0', '¦È=¦Ð/16','¦È=2¦Ð/16','¦È=3¦Ð/16','¦È=4¦Ð/16')
xlabel('percentage of kh')
ylabel('\epsilon (\theta)')
grid on