clear;clc;
close all
M=7;
h=20;
v=1500;

tau=0.001;
r=v*tau/h;
r=0;
%%%%%%%%%%%%%%%%%%%³õÊ¼ÏµÊý
c=zeros(1,M);
for m=1:M
    c(m)=1;
    temp=(-1)^(m+1)/(2*m-1);
    for n=1:M
        if n~=m
            c(m)=abs(((2*n-1)^2-r^2)/((2*n-1)^2-(2*m-1)^2))*c(m);
        end
    end
    c(m)=temp*c(m);
end



v=2598;%%% figure 1b
r=v*tau/h;

k=linspace((pi)/(100*h),(pi)/h,100);
for kk=1:5
    %     xita=(kkk-1)*15*pi/180;
    % xita=0
    % xita=pi/8
    %     xita=pi/4;
    xita=(kk-1)*pi/16;
    a=0;
    for n=1:M
        a=a+c(n)*sin((n-0.5)*k*h*cos(xita));
    end
    
    b=0;
    for n=1:M
        b=b+c(n)*sin((n-0.5)*k*h*sin(xita));
    end
    
    a=a.^2;
    b=b.^2;
    
    delta=2./(r*k*h);
    delta=delta.*asin(r*sqrt(a+b));
%     a1=(h/v*(1./delta-1));
    
    a1=delta;
    
    %     a1=a1*10^5;
    %     hold on
    %     plot(v*k*h/(2*pi*h),(a1),'k')
    
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
% axis([0 100 -1.5*10^-5  5*10^-5])
axis([0 100 0.97  1.01])
% axis([0 100 0.98  1.04])
legend('¦È=0', '¦È=¦Ð/16','¦È=2¦Ð/16','¦È=3¦Ð/16','¦È=4¦Ð/16')
xlabel('percentage of kh')
ylabel('\epsilon (\theta)')
grid on