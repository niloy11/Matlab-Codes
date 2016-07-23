A= 775;
sigbA = 3.7*10^9 ;
L= [ 5000 ]*10^-9;
H=0.35*10^-9;
%500  700 1000 3000 5000

km=.10;
Rb=10*10^-9;


f=0:0.001:0.1

for j=1:length(L)
for i=1:length(f)
    
    kp(j) =(sigbA )/(1/(L(j))+2/(pi *A));
    
    a11(j)=L(j);
    a21(j)=L(j);
    a31(j)=H;
    
    ak1=Rb*km
    alp1=Rb*km/H
    
    p(j)=H/L(j)

    L11(j)=((p(j))^2/(2*((p(j))^2-1)))+((p(j))/(2*(1-(p(j))^2)^(3/2)))*acos(p(j));
    gama(j)=(1+2*p(j))*alp1;

    L22(j)=L11(j);
    L33(j)=1-2*L11(j)

    kc11(j)=kp(j)/(1+gama(j)*L11(j)*kp(j)/km);
    kc22(j)=kp(j)/(1+gama(j)*L22(j)*kp(j)/km);
    kc33(j)=kp(j)/(1+gama(j)*L33(j)*kp(j)/km);

    b11(j)=(kc11(j)-km)/(km +L11(j)*(kc11(j)-km));
    b22(j)=(kc22(j)-km)/(km +L22(j)*(kc22(j)-km));
    b33(j)=(kc33(j)-km)/(km +L33(j)*(kc33(j)-km));

    
    
 
  k1(i) =(3+f(i)*(2*b11(j)*(1-L11(j))+b33(j)*(1-L33(j))))/(3-f(i)*(2*b11(j)*L11(j)+b33(j)*L33(j)));
  ka1(i)=k1(i)*km;
  ke1(i)=ka1(i)
  
  
end
f1=f*100;
plot(f1,ke1);
hold on;
xlabel('volume fraction f '), ylabel('thermal conductivity Kc (W/mK)'),
title('composite thermal conductivity of MLG-Silicone elastomer as variation with volume fraction for different lateral size ')
end
f11=f1';
ke11=ke1';