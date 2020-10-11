
% This script provides the ideal nozzle geometry using the method of
% of characteristics for an almost 2D diverging nozzle. Suppose gas is
% exhaust from a combustion chamber with no mass flow, and temperature and pressure are known.
% Using the 2D nozzle flow relationships, an optimal throat area is found
% that will produce the maximum amount of thrust for the given ambient pressure
% and parameters of the combustion chamber.
clear all
clc
T_c=2000;%Combustion chamber Temp
P_c=1.2*10^6;%Combustion chamber Pressure
P_amb=1.01*10^5;
T_amb=300;
gamma=1.25;
larghezza=0.4;
altezza=0.025/2;
R=287;
delta_err=0.00001;
delta_M=0.00001;
Tt=T_c/(1+(gamma-1)/2);
Pt=P_c/(1+(gamma-1)/2)^(gamma/(gamma-1));
rhot=Pt/R/Tt;
Vt=sqrt(gamma*R*Tt);

num=10;
theta_start=0.03;

Me=sqrt(((P_c/P_amb)^((gamma-1.)/gamma)-1)*2/(gamma-1))
Te=(P_amb/P_c)^((gamma-1)/gamma)*T_c;
rhoe=P_amb/R/Te;
Ve=Me*sqrt(gamma*R*Te);
% Ae=altezza*larghezza*(rhot*Vt)/(rhoe*Ve);
% Max_thrust=(P_amb/R/Te)*Ve*Ae*Ve;

f=(gamma+1)/(gamma-1);
nu_e=sqrt(f)*atan(sqrt(1/f*(Me^2-1)))-atan(sqrt(Me^2-1))
theta_max=(180/pi)*nu_e/2;
nu_A=theta_max;
nu_s=zeros(num);
theta_s=zeros(num);
for i=1:num
    nu_s(i)=((theta_max-theta_start)/(num-1))*(i-1)+theta_start;
    theta_s(i)=nu_s(i);
end

theta=zeros(num);
nu=zeros(num);

nu(1,1)=2*theta_start;
theta(1,1)=0.0;
for i=2:num
    for j=1:i
        
        if (j==1)
            theta(i,j)=[theta_s(i)+nu_s(i)  +  theta(i-1,j)-nu(i-1,j)]/2;
            nu(i,j)=[theta_s(i)+nu_s(i)  -  (theta(i-1,j)-nu(i-1,j))]/2;
        
        elseif (j==i)
            theta(i,j)=0;
            nu(i,j)=theta(i,j-1)+nu(i,j-1);    
        
        elseif(j~=i)
             theta(i,j)=[theta(i,j-1)+nu(i,j-1)  +  theta(i-1,j)-nu(i-1,j)]/2;
             nu(i,j)=[theta(i,j-1)+nu(i,j-1)  -  theta(i-1,j)+nu(i-1,j)]/2;
           
       
        end
    end
end


nu_2=0.0;
error=100;
MM_1=0.0;
MM=0.0;
mach=zeros(num);
mi=zeros(num);
error=100;
nu_1=0.0;
aa=sqrt((gamma+1)/(gamma-1));
bb=(gamma-1)/(gamma+1);
mi_s=zeros(num)
for i=1:num
      nu_is=(pi/180)*nu_s(i);
    
while(error>0.0001)
        nu_2=aa*atan(sqrt(bb*(MM_1^2-1)))-atan(sqrt(MM_1^2-1));
        error=nu_is-nu_2;
        MM_1=MM_1+0.0001;
end
mi_s(i)=(180/pi)*asin(1/MM_1);

error=100;
for j=1:i
   nu_ii=(pi/180)*nu(i,j);
   while(error>0.0001)
        nu_1=aa*atan(sqrt(bb*(MM^2-1)))-atan(sqrt(MM^2-1));
        error=nu_ii-nu_1;
        MM=MM+0.0001;
end
mach(i,j)=MM;
mi(i,j)=(180/pi)*asin(1/MM);
MM=0.0;
MM_1=0.0;
error=100;
end
end


angular_plus=theta+mi;
angular_minus=theta-mi;
theta_m=theta_max;
P(1,1,:)=intersez(tan((pi/180)*((theta_s(1)-mi_s(1)))),[0 altezza],tan((pi/180)*(theta_s(1)+mi_s(1))),[0 -altezza]);
  
  for i=2:num
      P(i,1,:)=intersez(tan((pi/180)*((theta_s(i)-mi_s(i)))),[0 altezza],tan(pi/180)*(angular_plus(i-1,1)+2*theta_max),P(i-1,1,:));
  end
  
  for car=2:num
         P(car,car,:)=intersez(tan((pi/180)*angular_minus(car,car-1)),P(car,car-1,:),0,[0 0]);
         for i=(car+1):num
              P(i,car,:)=intersez(tan((pi/180)*angular_minus(i,car-1)),P(i,car-1,:),tan((pi/180)*angular_plus(i-1,car)),P(i-1,car,:));
         end
  end

    P_c(1,1)=0;
    P_c(1,2)=altezza;
    m=theta_max;
    
    for i=2:num+1
        P_c(i,:)=intersez(tan((pi/180)*m),P_c(i-1,:),tan((pi/180)*angular_plus(num,i-1)),P(num,i-1,:));
        m=theta(num,i-1);
    end
    
    
    A1=altezza*2*larghezza
    A2=P_c(i,2)*2*larghezza
    rapp=A1/A2
    
    
    
for i=1:num
    xx(i)=P(i,i,1);
    yy(i)=P(i,i,2);
end   
P(:,:,1)=P(:,:,1)+transpose(P(:,:,1));
P(:,:,2)=P(:,:,2)+transpose(-P(:,:,2));
for i=1:num
    P(i,i,1)=xx(i);
    P(i,i,2)=yy(i);
end
    l=1;
    k=1;
     for i=1:num
         for j=1:num
             PP(j,i,1)=P(j,i,1);
             PP(j,i,2)=P(j,i,2);
             
             if j==num
             PP(j+1,i,1)=P_c(i+1,1);
             PP(j+1,i,2)=P_c(i+1,2);
             end
          end
             
        
     end
    for i=1:num
        
        fan(i,1,1)=PP(i,1,1);
        fan(i,1,2)=PP(i,1,2);
        fan(i,2,1)=0;
        fan(i,2,2)=altezza;
    end
    
    figure(3)
    hold on
    plot(0,altezza,'*')
    title('Supersonic nozzle')
    xlabel('lenth [m]')
    ylabel('width [m]')
    legend(sprintf('n characteristics: %d\nExit mach value:%f',num,Me))
    plot(transpose(fan(:,:,1)),transpose(fan(:,:,2)),'red')
    plot(PP(:,:,1),PP(:,:,2),'red')
    plot(P_c(:,1),P_c(:,2),'*')
    plot(P_c(:,1),P_c(:,2))
    plot(0,-altezza,'*')
    plot(transpose(fan(:,:,1)),-transpose(fan(:,:,2)),'red')
    plot(PP(:,:,1),-PP(:,:,2),'red')
    plot(P_c(:,1),-P_c(:,2),'*')
    plot(P_c(:,1),-P_c(:,2))
    hold off
    A1=altezza*2*larghezza
    A2=P_c(i,2)*2*larghezza
    rapp=A2/A1