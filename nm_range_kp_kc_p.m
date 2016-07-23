clc;
close all;
clear all;


%%%%%%%%%%%
A= 775;
sigbA = 3.7*10^9 ;
% p=[ 0.001167 0.00233 .0035 .005 0.01 ]
p=[ 0.0098  ]
%0.0007 0.0021 .0049  0.0098

L= 500*10^-9

% v=2*10^4;
% C= 700 *2300

%kp= 2000;
km=0.1;

Rb=10*10^-9;


f=0:0.001:0.1

for j=1:length(p)
for i=1:length(f)
    
    Lef=(A* L)/(A+ L);
    
%     kp = 1/3 *Lef* v* C*0.5
    kp =(sigbA )/(1/(L)+2/(pi *A));
    
    H(j)=p(j)*L
    
    a11(j)=L;
    a21(j)=L;
    a31(j)=H(j);
    
    ak1=Rb*km
    alp1(j)=Rb*km/H(j)
    
    

    L11(j)=((p(j))^2/(2*((p(j))^2-1)))+((p(j))/(2*(1-(p(j))^2)^(3/2)))*acos(p(j))
    gama(j)=(1+2*p(j))*alp1(j);

    L22(j)=L11(j);
    L33(j)=1-2*L11(j)

    kc11(j)=kp/(1+gama(j)*L11(j)*kp/km);
    kc22(j)=kp/(1+gama(j)*L22(j)*kp/km);
    kc33(j)=kp/(1+gama(j)*L33(j)*kp/km);

    b11(j)=(kc11(j)-km)/(km +L11(j)*(kc11(j)-km));
    b22(j)=(kc22(j)-km)/(km +L22(j)*(kc22(j)-km));
    b33(j)=(kc33(j)-km)/(km +L33(j)*(kc33(j)-km));

    
    
 
  k1(i) =(3+f(i)*(2*b11(j)*(1-L11(j))+b33(j)*(1-L33(j))))/(3-f(i)*(2*b11(j)*L11(j)+b33(j)*L33(j)));
  ka1(i)=k1(i)*km;
  ke1(i)=ka1(i);
  
  
end
plot(f*100,ke1)
hold on;
xlabel('Volume fraction f '), ylabel('Thermal conductivity Kc (W/mK)'),
title('composite thermal conductivity of MLG-Silicone elastomer vs variation with volume fraction for different a ')
end
f1=f*100;
f2=f1';
ke2=ke1';