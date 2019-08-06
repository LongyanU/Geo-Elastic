function F = myfun2(x)
global M v h dt




k=linspace(1/1000000,0.82*pi/h,25);



r=v*dt/h;
F=zeros(5,25);

for jj=1:5
    
    xita=(jj-1)*pi/16;
    
    temp1=0;
    for m=1:M
        temp1=temp1+x(m)*sin((m-1/2)*k*sin(xita)*h);
    end
    temp=temp1.^2;
    
    temp2=0;
    for m=1:M
        temp2=temp2+x(m)*sin((m-1/2)*k*cos(xita)*h);
    end
    temp=temp+temp2.^2;
    
    bb=1/(2*r^2)*(-cos(k*r*h)+1);
    F(jj,:)=temp./bb;
end

F=F-1;
F=F(1,:)*F(1,:)'+F(2,:)*F(2,:)'+F(3,:)*F(3,:)'+F(4,:)*F(4,:)'+F(5,:)*F(5,:)';