clc
clear all 
close all
L= [1 3 5]*10^-6;
H3=0.35*10^-9;
vt=18000;
vl=24000;
km=0.1
f=0:0.001:0.1

for j= 1:length(L)
for i=1:length(f)
    
gammat=0.75 ;
gammal= 1.8;
M=1.87*1.992*10^-26;
wmaxt=108*10^12;
wmaxl=241*10^12;
kb=1.38065*10^-23;

hc=1.0545*10^-34;
h=0.35*10^-9;

T=300;

        wmint(j)=vt/gammat*sqrt(M*vt*wmaxt/(kb*T*L(j)));
        wminl(j)=vl/gammal*sqrt(M*vl*wmaxl/(kb*T*L(j)));
        Ft(j)=-log(abs(exp(hc*wmint(j)/(kb*T))-1))+hc*wmint(j)/(kb*T)*(exp(hc*wmint(j)/(kb*T))/(exp(hc*wmint(j)/(kb*T))-1));
        Fl(j)=-log(abs(exp(hc*wminl(j)/(kb*T))-1))+hc*wminl(j)/(kb*T)*(exp(hc*wminl(j)/(kb*T))/(exp(hc*wminl(j)/(kb*T))-1));
        kp(j)=M/(4*pi*T*h)*(wmaxt*vt^2/((gammat)^2)*Ft(j)+wmaxl*vl^2/((gammal)^2)*Fl(j));
        
    %%%%%

 a1(j)=L(j);
 a2(j)=L(j);
 a3(j)=H3;


Rb=10*10^-9;

p(j)=H3/L(j);

ak3=Rb*km;
alp3=Rb*km/H3;


L11(j)=((p(j))^2/(2*((p(j))^2-1)))+((p(j))/(2*(1-(p(j))^2)^(3/2)))*acos(p(j));
gama(j)=(1+2*p(j))*alp3;

L22(j)=L11(j);
L33(j)=1-2*L11(j);

kc11(j)=kp(j)/(1+gama(j)*L11(j)*kp(j)/km);
kc22(j)=kp(j)/(1+gama(j)*L22(j)*kp(j)/km);
kc33(j)=kp(j)/(1+gama(j)*L33(j)*kp(j)/km);

% kc113=kp3 +gama3*(1-L113)*km;
% kc223=kp3 +gama3*(1-L223)*km;
% kc333=kp3 +gama3*(1-L333)*km;

b11(j)=(kc11(j)-km)/(km +L11(j)*(kc11(j)-km));
b22(j)=(kc22(j)-km)/(km +L22(j)*(kc22(j)-km));
b33(j)=(kc33(j)-km)/(km +L33(j)*(kc33(j)-km));  
    %%%%%    

  kc1(i) =(3+f(i)*(2*b11(j)*(1-L11(j))+b33(j)*(1-L33(j))))/(3-f(i)*(2*b11(j)*L11(j)+b33(j)*L33(j)));
  ktim(i)=kc1(i)*km;
  
end

 plot(f*100,ktim)
 hold on;
xlabel('Volume fraction'), ylabel('thermal conductivity Kc (W/mK)'),
% title('composite thermal conductivity of MLG-Silicone elastomer as variation with volume fraction for different lateral size ')
end
f1=(f*100);
f2=f1';
ktim1=ktim';