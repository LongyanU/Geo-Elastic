function F = myfun3(x)

global M  ratio h
k=linspace(1/1000000,ratio*pi/h,100);
temp1=0;
for m=1:M
    temp1=temp1+2*x(m)*sin((m-1/2)*k*h);
end

temp=temp1.^2;

F=(temp-k.^2*h^2)*(temp-k.^2*h^2)';