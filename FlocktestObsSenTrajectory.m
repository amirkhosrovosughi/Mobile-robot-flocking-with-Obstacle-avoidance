%%Amirkhosro Vosughi 11463709 WSU
%Research Project for course EE505
%Simulstion of "Distributed Flocking Control of Mobile Robots by Bounded Feedback"

clear all
clc

%% Initialization
n1 = 3; %==>> number of robots
n2 = 2; %==>> number of onstab
n=n1+n2;
% 
 
Group = [ones(1,n1),2*ones(1,n2)]';  %==>> Determines with robot determines to each group

%Robot parameters
d= 0.1;
r=5; %==>> communication range
r1=1; %==>> obstacle sensing range

%controller parameters
Vref=1;
Wref=0.25;

b11=0;
b21=0;
b31=0;

b12=0;
b22=0;
b32=0;

%Generate initial location
%obstable
 xObs=[];
 yObs=[];

xObs=[2,-4]';
yObs=[-4,-10]';
% xObs=[2,-1]';
% yObs=[1,-2]';
%robots inital condition
x0= [0,1.5,4.5]';
y0= [0,2,0]';
   
%total vecotrs
if (n2~=0)
 x0=[x0;xObs];
 y0=[y0;yObs];
end

teta_init=zeros(n1,1);
for k=1:length(teta_init)
    teta_init(k)=(rand(1));
end

teta_init =[0.6128
    0.8194
    0.5319];

teta_init = [teta_init;zeros(n2,1)];
%teta_init=[0,0]';

v_init=zeros(n1,1);
for k=1:length(v_init)
    v_init(k)=1*(rand(1)+1);
end

v_init =[1.2021
    1.4539
    1.4279];

v_init = [v_init;zeros(n2,1)];
%v_init=[1,1]';

w_init=zeros(n1,1);
for k=1:length(w_init)
    w_init(k)=1*(rand(1)+1);
end

w_init =[1.9661
    1.6201
    1.6954];
w_init = [w_init;zeros(n2,1)];
%w_init=[0.2,0.2]';
%% Simulation
T=0.01; %sec
Tfinal=400; %sec
T_vec=0:T:Tfinal;
Tn=length(T_vec);

%%Definition of process variables

u1= zeros(n,Tn);
u2= zeros(n,Tn); %===>> select size n not n1, because maybe we want to have moving obstable in future

v=zeros(n,Tn);
v(:,1)=v_init;    %%===>> select size n not n1, because maybe we want to have moving obstable in future

w=zeros(n,Tn);
w(:,1)=w_init; %===>> select size n not n1, because maybe we want to have moving obstable in future

x=zeros(n,Tn);
x(:,1)=x0; %===>> select size n not n1, because maybe we want to have moving obstable in future

y=zeros(n,Tn);
y(:,1)=y0; %===>> select size n not n1, because maybe we want to have moving obstable in future

teta=zeros(n,Tn);
teta(:,1)=teta_init;  %===>> select size n not n1, because maybe we want to have moving obstable in future


N=zeros(n,n); 
A=zeros(n,n);
%L=zeros(n,n);

xdot=zeros(n,1);
ydot=zeros(n,1);

tempoli=zeros(n,Tn);
tempoli2=zeros(n,Tn);

%% Simulation
for k=2:Tn

%Update The graphs
N=zeros(n,n); 
A=zeros(n,n);

for i=1:n
    for j=i:n
        temp1=x(i,k-1)-x(j,k-1);
        temp2=y(i,k-1)-y(j,k-1);
        temp3=(temp1^2+temp2^2)^0.5;
        
        if (temp3<r && i~=j)
            N(i,j)=1;
            N(j,i)=1;
            
            A(i,j)=temp3;
            A(j,i)=temp3;
        end     
    end
end

%Calculate control input
% calculate u(k)

%calculate u1(k) ==>> assumption of understanding of gradient V to pij
for i=1:n1
    
    sigma1=0;
    for j=1:n
        if (N(i,j)==1 && Group(i)==Group(j) )  
            temp1=A(i,j)*(v(i,k-1)-v(j,k-1));
            sigma1=sigma1+temp1;
        end
    end
    
    sigma2=0;
    for j=1:n
        if (N(i,j)==1 && Group(i)==Group(j) ) 
            temp1=xdot(i)-xdot(j);
            temp2=ydot(i)-ydot(j);
            temp3=A(i,j)*(temp1*cos(teta(i,k-1))+temp2*sin(teta(i,k-1)));
            sigma2= sigma2+temp3;    
        end
    end
    %sigma2=0;
    
    sigma3=0;
    for j=1:n
        if (N(i,j)==1) 
            
            temp1=0;
            temp4=0;
            if ( Group(i)==Group(j))
                temp1=-2*A(i,j)^-3+2*A(i,j)/(r^2-A(i,j)^2)+2*b11*A(i,j)+b21;
                
            elseif (A(i,j)<r1)
                temp1=(-2*A(i,j)^-3+2*b12*A(i,j)+b22);
               
            end
            temp2= temp1*((x(i,k-1)-x(j,k-1))/A(i,j))*cos(teta(i,k-1))+temp1*((y(i,k-1)-y(j,k-1))/A(i,j))*sin(teta(i,k-1));
              
            tempoli(i,k)= temp2;
            
            
            sigma3= sigma3+ temp2;    
        end
    end
    %sigma3=0;
    
    u1(i,k)= -sigma1-sigma2-sigma3;
end

%%%%%%%%%%%%%%trajectoy for leader robot
temp1=Vref-v(1,k-1);
u1(1,k)=u1(1,k)+temp1*0.2;
%%%%%%%%%%%%%%%
% calculate u2(k)

for i=1:n1
    
    sigma1=0;
    for j=1:n
        if (N(i,j)==1 && Group(i)==Group(j) )  
            temp1=A(i,j)*(w(i,k-1)-w(j,k-1));
            sigma1=sigma1+temp1;
        end
    end
    
    sigma2=0;
    for j=1:n
        if (N(i,j)==1 && Group(i)==Group(j) ) 
            temp1=xdot(i)-xdot(j);
            temp2=ydot(i)-ydot(j);
            temp3=A(i,j)*(-temp1*sin(teta(i,k-1))+temp2*cos(teta(i,k-1)));
            sigma2= sigma2+temp3;    
        end
    end
    
    sigma3=0;
    for j=1:n
        if (N(i,j)==1) 
            
            temp1=0;
            if ( Group(i)==Group(j))
                temp1=-2*A(i,j)^-3+2*A(i,j)/(r^2-A(i,j)^2)+2*b11*A(i,j)+b21;
            elseif (A(i,j)<r1)
                temp1=-2*A(i,j)^-3+2*b12*A(i,j)+b22;
            end
            temp2= -temp1*((x(i,k-1)-x(j,k-1))/A(i,j))*sin(teta(i,k-1))+temp1*((y(i,k-1)-y(j,k-1))/A(i,j))*cos(teta(i,k-1));
            
            sigma3= sigma3+temp2;    
        end
    end
    %sigma3=0;
    
    u2(i,k)= -sigma1-sigma2-sigma3;
end

%%%%%%%%%%%%%%trajectoy for leader robot
temp1=Wref-w(1,k-1);
u2(1,k)=u2(1,k)+temp1*0.2;
%%%%%%%%%%%%%%%
%% Dynamic of system


for j=1:n
    v(j,k)=v(j,k-1)+u1(j,k)*T;
    w(j,k)=w(j,k-1)+(u2(j,k)/d)*T;
end
for j=1:n
   teta(j,k)=teta(j,k-1)+w(j,k)*T;
   
   xdot(j)=v(j,k)*cos(teta(j,k-1))-d*w(j,k)*sin(teta(j,k-1));
   
   ydot(j)=v(j,k)*sin(teta(j,k-1))+d*w(j,k)*cos(teta(j,k-1));

   x(j,k)=x(j,k-1)+xdot(j)*T;
   y(j,k)=y(j,k-1)+ydot(j)*T;
end
    

end

%% Plotting the results

figure(1)
plot(x(1,:),y(1,:))
hold on

plot(x(2,:),y(2,:))
hold on


plot(x(3,:),y(3,:))
hold on



% plot(x(4,:),y(4,:))
% hold on

plot(xObs,yObs,'r*')

hold off

% plot(x(5,:),y(5,:))
% hold off



% figure(2)
% plot(x(6,:),y(6,:))
% hold on
% 
% plot(x(7,:),y(7,:))
% hold on
% 
% plot(x(8,:),y(8,:))
% hold on
% 
% plot(x(9,:),y(9,:))
% hold on
% 
% plot(x(10,:),y(10,:))
% hold off