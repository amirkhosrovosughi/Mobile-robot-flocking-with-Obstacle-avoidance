clear all
r=5;
 x=0.4:0.0001:4.99; %==>> good to display
%x=0.2:0.01:6.9;
n=length(x);

%e=0.05

for i=1:n
fx(i)=-2*x(i)^-3;
gx(i)=2*x(i)/(r^2-x(i)^2);
Ax(i)=fx(i)+gx(i);

Bx(i)=x(i)^-2-log(r^2-x(i)^2);
Cx(i)=x(i)^-2-0.0400;

%Axrx=Ax(i)/x(i);
end


plot(x,Cx)
xlim([0.4,4.99])