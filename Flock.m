%%Amirkhosro Vosughi 11463709 WSU
%Research Project for course EE505
%Simulstion of "Distributed Flocking Control of Mobile Robots by Bounded Feedback"

clear all
clc

%% Initialization
n = 10; %==>> To generelize code later
% robots 1 to 5 belongs to group 1
% robots 6 to 10 belongs to group 2
Group = [1,1,1,1,1,2,2,2,2,2]';  %==>> Determines with robot determines to each group
%Group = [1,1,2,1,1,2,2,2,1,2]';

%Robot parameters
d= 0.1;
r=5;

%controller parameters
b11=0;
b21=0;
b31=0;

b12=0;
b22=0;
b32=0;

%Generate initial location randomly
% x0= [0,-2,2,-2,-1,1,1.5,2,-1.5,1]';
% y0= [0,2,2,5,6,6,6,6,8,8]';
    
x0= [-2,-2,-1.5,-1,0,1,1,1.5,2,2]';
y0= [5,2,8,6,0,8,6,6,6,2]';


teta_init=zeros(n,1);
for k=1:length(teta_init)
    teta_init(k)=(rand(1));
end

teta_init=[0.7923, 0.6988, 0.5291, 0.7695, 0.1557, 0.8052, 0.7741, 0.8293, 0.7723, 0.5643]';

v_init=zeros(n,1);
for k=1:length(v_init)
    v_init(k)=1*(rand(1)+1);
end

v_init=[0.6222, 1.8468, 0.8604, 0.3696, 1.8098, 1.9595, 0.8777, 0.2222, 0.5161, 0.8174]';

w_init=zeros(n,1);
for k=1:length(w_init)
    w_init(k)=1*(rand(1)+1);
end

w_init=[0.2974, 0.1311, 0.3014, 0.3556, 0.1109, 0.0587, 0.1483, 0.1594, 0.2121, 0.2539]';

%% Simulation
T=0.01; %sec
Tfinal=10; %sec
T_vec=0:T:Tfinal;
Tn=length(T_vec);

%%Definition of process variables

u1= zeros(n,Tn);
u2= zeros(n,Tn);

v=zeros(n,Tn);
v(:,1)=v_init;

w=zeros(n,Tn);
w(:,1)=w_init;

x=zeros(n,Tn);
x(:,1)=x0;

y=zeros(n,Tn);
y(:,1)=y0;

teta=zeros(n,Tn);
teta(:,1)=teta_init;


N=zeros(n,n); 
A=zeros(n,n);
L=zeros(n,n);


xdot=zeros(n,1);
ydot=zeros(n,1);

tempoli=zeros(n,Tn);
tempoli2=zeros(n,Tn);

minA=zeros(1,Tn);
maxA=zeros(1,Tn);
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
  
minA(k)=min(min(A));
maxA(k)=max(max(A));
%Calculate control input
% calculate u(k)

%calculate u1(k) ==>> assumption of understanding of gradient V to pij
for i=1:n
    
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
                %temp4=-A(i,j)^-2+1/(r^2-A(i,j));%+2*b11*A(i,j)+b21;
                temp4=-A(i,j)^-4+1/(r^2-A(i,j)^2);%+2*b11*A(i,j)+b21;
            else
                temp1=-2*A(i,j)^-3+2*b12*A(i,j)+b22;
                %temp4=-A(i,j)^-2;
                temp4=-A(i,j)^-4;
            end
            temp2= temp1*((x(i,k-1)-x(j,k-1))/A(i,j))*cos(teta(i,k-1))+temp1*((y(i,k-1)-y(j,k-1))/A(i,j))*sin(teta(i,k-1));
            temp3= temp4*((x(i,k-1)-x(j,k-1))*2)*cos(teta(i,k-1))+temp4*((y(i,k-1)-y(j,k-1))*2)*sin(teta(i,k-1));
            
            %temp2= temp1;
            tempoli(i,k)= temp2;
            tempoli2(i,k)= temp3;
            
            sigma3= sigma3+ temp2;    
        end
    end
    %sigma3=0;
    
    u1(i,k)= -sigma1-sigma2-sigma3;
end


%%%%%%%%%%%%%%%%%%not edited from here %%%%%%%%%%

% calculate u2(k)

for i=1:n
    
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
            else
                temp1=-2*A(i,j)^-3+2*b12*A(i,j)+b22;
            end
            temp2= -temp1*((x(i,k-1)-x(j,k-1))/A(i,j))*sin(teta(i,k-1))+temp1*((y(i,k-1)-y(j,k-1))/A(i,j))*cos(teta(i,k-1));
            
            sigma3= sigma3+temp2;    
        end
    end
    %sigma3=0;
    
    u2(i,k)= -sigma1-sigma2-sigma3;
end

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

plot(x(4,:),y(4,:))
hold on

plot(x(5,:),y(5,:))
%hold off



%figure(2)
plot(x(6,:),y(6,:))
hold on

plot(x(7,:),y(7,:))
hold on

plot(x(8,:),y(8,:))
hold on

plot(x(9,:),y(9,:))
hold on

plot(x(10,:),y(10,:))
hold off