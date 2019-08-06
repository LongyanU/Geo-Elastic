clear;
clc;
close all
global M v h dt

dt=0.001;
h=20;
v=1500;

options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);
lb=-1.5*ones(M,1);
ub=1.5*ones(M,1);
tic
M=7;
x0=0.01*ones(1,M);


[x,fval,out,iteration]= fmincon(@myfun2,x0,[],[],[],[],lb,ub);    % Invoke optimizer
%if it does not work, please run the code again. fmincon do not make get
%the optimized FD coefficient every time.

v=2598 ; %% figure 2b
r=v*dt/h;
c=x
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
    
    a=a.*a;
    b=b.*b;
    
    delta=2./(r*k*h);
    delta=delta.*asin(r*sqrt(a+b));

    
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
