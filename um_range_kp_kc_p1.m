

clc;
clear all;
close all;

L= 5*10^-6;

H=[ 4 ]*0.35 *10^-9;
% H3=3*0.35*10^-9;
vt=18000;
vl=24000;
km=0.1

 kp= [ 4000 2700 2200 1300]
% kp= [ 2200 2200 2200 2200]

f=0:0.001:0.1

for j= 1:length(H)
for i=1:length(f)
    
% gammat=0.75 ;
% gammal= 1.8;
% M=1.87*1.992*10^-26;
% wmaxt=108*10^12;
% wmaxl=241*10^12;
% kb=1.38065*10^-23;
% 
% hc=1.0545*10^-34;
% h(j)=p(j)*L;
% % T=250:1:500;
% T=300;
% 
%         wmint=vt/gammat*sqrt(M*vt*wmaxt/(kb*T*L));
%         wminl=vl/gammal*sqrt(M*vl*wmaxl/(kb*T*L));
%         Ft=-log(abs(exp(hc*wmint/(kb*T))-1))+hc*wmint/(kb*T)*(exp(hc*wmint/(kb*T))/(exp(hc*wmint/(kb*T))-1));
%         Fl=-log(abs(exp(hc*wminl/(kb*T))-1))+hc*wminl/(kb*T)*(exp(hc*wminl/(kb*T))/(exp(hc*wminl/(kb*T))-1));
%         kp(j)=M/(4*pi*T*h(j))*(wmaxt*vt^2/((gammat)^2)*Ft+wmaxl*vl^2/((gammal)^2)*Fl);
        
    %%%%%
 
 a1(j)=L;
 a2(j)=L;
 a3(j)=H(j);


Rb=10*10^-9;

p(j)=H(j)/L;

% ak3=Rb*km;
alp3(j)=Rb*km/H(j);


L11(j)=((p(j))^2/(2*((p(j))^2-1)))+((p(j))/(2*(1-(p(j))^2)^(3/2)))*acos(p(j));
gama(j)=(1+2*p(j))*alp3(j);

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
xlabel('volume fraction f '), ylabel('thermal conductivity Kc (W/mK)'),
% title('composite thermal conductivity of MLG-Silicone elastomer vs variation with volume fraction for different a ')

end
f2=(f*100)';
ktim1=ktim';
% figure(2)
% plot(L,kp)

%%%%%%%%
%%%%%%%%

% figure(3)
% 
% T=300:400
%  f=0.1
% 
% for j= 1:length(L)
% for i=1:length(T)
%     
% gammat=0.75 ;
% gammal= 1.8;
% M=1.87*1.992*10^-26;
% wmaxt=108*10^12;
% wmaxl=241*10^12;
% kb=1.38065*10^-23;
% 
% hc=1.0545*10^-34;
% h=3*0.35*10^-9;
% 
% 
%         wmint(i)=vt/gammat*sqrt(M*vt*wmaxt/(kb*T(i)*L(j)));
%         wminl(i)=vl/gammal*sqrt(M*vl*wmaxl/(kb*T(i)*L(j)));
%         Ft(i)=-log(abs(exp(hc*wmint(i)/(kb*T(i)))-1))+hc*wmint(i)/(kb*T(i))*(exp(hc*wmint(i)/(kb*T(i)))/(exp(hc*wmint(i)/(kb*T(i)))-1));
%         Fl(i)=-log(abs(exp(hc*wminl(i)/(kb*T(i)))-1))+hc*wminl(i)/(kb*T(i))*(exp(hc*wminl(i)/(kb*T(i)))/(exp(hc*wminl(i)/(kb*T(i)))-1));
%         kp(i)=M/(4*pi*T(i)*h)*(wmaxt*vt^2/((gammat)^2)*Ft(i)+wmaxl*vl^2/((gammal)^2)*Fl(i));
%         
%     %%%%%
% 
%  a1(j)=L(j);
%  a2(j)=L(j);
%  a3(j)=H3;
% 
% 
% Rb=5*10^-9;
% 
% p(j)=H3/L(j);
% 
% ak3=Rb*km;
% alp3=Rb*km/H3;
% 
% 
% L11(j)=((p(j))^2/(2*((p(j))^2-1)))+((p(j))/(2*(1-(p(j))^2)^(3/2)))*acos(p(j));
% gama(j)=(1+2*p(j))*alp3;
% 
% L22(j)=L11(j);
% L33(j)=1-2*L11(j);
% 
% kc11(i)=kp(i)/(1+gama(j)*L11(j)*kp(i)/km);
% kc22(i)=kp(i)/(1+gama(j)*L22(j)*kp(i)/km);
% kc33(i)=kp(i)/(1+gama(j)*L33(j)*kp(i)/km);
% 
% % kc113=kp3 +gama3*(1-L113)*km;
% % kc223=kp3 +gama3*(1-L223)*km;
% % kc333=kp3 +gama3*(1-L333)*km;
% 
% b11(i)=(kc11(i)-km)/(km +L11(j)*(kc11(i)-km));
% b22(i)=(kc22(i)-km)/(km +L22(j)*(kc22(i)-km));
% b33(i)=(kc33(i)-km)/(km +L33(j)*(kc33(i)-km));  
%     %%%%%    
% 
%   kc1(i) =(3+f*(2*b11(i)*(1-L11(j))+b33(i)*(1-L33(j))))/(3-f*(2*b11(i)*L11(j)+b33(i)*L33(j)));
%   ktim2(i)=kc1(i)*km;
%   
% end
% 
%  plot(T,ktim2)
%  hold on;
% 
% end


