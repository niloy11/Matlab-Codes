clc;
close all;
clear all;


%%%%%%%%%%%
A1= 500 *10^-9;
sigbA = 3.7*10^9 ;

L1= 300*10^-9;%%%%%%%%%%

% v=2*10^4;
% C= 700*2300  
%  kp1 = 1/3 *Lef1* v* C
 
 kp1 =(sigbA )/(1/(L1)+2/(pi *A1));
%kp= 2000;
km=0.1;
H1=3*.35*10^-9;%%%%%%%%%%%
% H=.35*10^-9
a11=L1;
a21=L1;
a31=H1;

% kp1=550%%%%%%%%%%

Rb=10*10^-9;

p1=H1/L1

ak1=Rb*km
alp1=Rb*km/H1


L111=(p1^2/(2*(p1^2-1)))+(p1/(2*(1-p1^2)^(3/2)))*acos(p1)
gama1=(1+2*p1)*alp1;

L221=L111;
L331=1-2*L111

kc111=kp1/(1+gama1*L111*kp1/km);
kc221=kp1/(1+gama1*L221*kp1/km);
kc331=kp1/(1+gama1*L331*kp1/km);

% kc111=kp1 +gama1*(1-L111)*km;
% kc221=kp1 +gama1*(1-L221)*km;
% kc331=kp1 +gama1*(1-L331)*km;

b111=(kc111-km)/(km +L111*(kc111-km));
b221=(kc221-km)/(km +L221*(kc221-km));
b331=(kc331-km)/(km +L331*(kc331-km));


%%%%%%%%%%%
%%%%%%%%%%%
A2=775*10^-9;
L2=100*10^-9;
H2=0.35*10^-9;

% H=.35*10^-9
a12=L2;
a22=L2;
a32=H2;

kp2 =(sigbA )/(1/(L2)+2/(pi *A2));

p2=H2/L2


ak2=Rb*km
alp2=Rb*km/H2

L112=(p2^2/(2*(p2^2-1)))+(p2/(2*(1-p2^2)^(3/2)))*acos(p2)
gama2=(1+2*p2)*alp2;

L222=L112;
L332=1-2*L112

kc112=kp2/(1+gama2*L112*kp2/km);
kc222=kp2/(1+gama2*L222*kp2/km);
kc332=kp2/(1+gama2*L332*kp2/km);

% kc112=kp2 +gama2*(1-L112)*km;
% kc222=kp2 +gama2*(1-L222)*km;
% kc332=kp2 +gama2*(1-L332)*km;

b112=(kc112-km)/(km +L112*(kc112-km));
b222=(kc222-km)/(km +L222*(kc222-km));
b332=(kc332-km)/(km +L332*(kc332-km));




%%%%%%%%%%%
L3= 5000*10^-9;%%%%%%%%%%
H3=3*0.35*10^-9;
% Lef3=(A* L3)/(A+ L3);
% v=1*10^6;
% C= 700
   kp3=2200
   
%    T=300
%     vt=18000;
%     vl=24000;
%     gammat=0.75;
%     gammal=1.8;
%     M=1.87*1.992*10^-26;
%     wmaxt=108*10^12;
%     wmaxl=241*10^12;
%     kb=1.38065*10^-23;
%     hc=1.0545*10^-34;
%     h=0.35*10^-9;
%         wmint=vt/gammat*sqrt(M*vt*wmaxt/(kb*T*L3));
%         wminl=vl/gammal*sqrt(M*vl*wmaxl/(kb*T*L3));
%         Ft=-log(abs(exp(hc*wmint/(kb*T))-1))+hc*wmint/(kb*T)*(exp(hc*wmint/(kb*T))/(exp(hc*wmint/(kb*T))-1));
%         Fl=-log(abs(exp(hc*wminl/(kb*T))-1))+hc*wminl/(kb*T)*(exp(hc*wminl/(kb*T))/(exp(hc*wminl/(kb*T))-1));
%         kumeq=M/(4*pi*T*h)*(wmaxt*vt^2/gammat^2*Ft(i)+wmaxl*vl^2/gammal^2*Fl(i));
%  
%         kp3(i)=kumeq(i)*1     ; 
   
   

a13=L3;
a23=L3;
a33=H3;

% kp3=2000%%%%%%%%%%%%%

p3=H3/L3

ak3=Rb*km
alp3=Rb*km/H3


L113=(p3^2/(2*(p3^2-1)))+(p3/(2*(1-p3^2)^(3/2)))*acos(p3)
gama3=(1+2*p3)*alp3;

L223=L113;
L333=1-2*L113

kc113=kp3/(1+gama3*L113*kp3/km);
kc223=kp3/(1+gama3*L223*kp3/km);
kc333=kp3/(1+gama3*L333*kp3/km);

% kc113=kp3 +gama3*(1-L113)*km;
% kc223=kp3 +gama3*(1-L223)*km;
% kc333=kp3 +gama3*(1-L333)*km;

b113=(kc113-km)/(km +L113*(kc113-km));
b223=(kc223-km)/(km +L223*(kc223-km));
b333=(kc333-km)/(km +L333*(kc333-km));



fr1=0.75
fr2=0.1
fr3=0.15





f=0:0.001:0.1
% cs=1/3; %%%%
for i=1:length(f)
%   k11(i)=(2+ f(i)*(b11*(1-L11)*(1+cs)+ b33*(1-L33)*(1-cs) ))/(2-f(i)*(b11*L11*(1+cs)+b33*L33*(1-cs)));
 k1(i) =(3+f(i)*(2*b111*(1-L111)+b331*(1-L331)))/(3-f(i)*(2*b111*L111+b331*L331));
  ka1(i)=k1(i)*km;
  ke1(i)=ka1(i)*fr1;
  k2(i) =(3+f(i)*(2*b112*(1-L112)+b332*(1-L332)))/(3-f(i)*(2*b112*L112+b332*L332));
  ka2(i)=k2(i)*km;
  ke2(i)=ka2(i)*fr2;
   k3(i) =(3+f(i)*(2*b113*(1-L113)+b333*(1-L333)))/(3-f(i)*(2*b113*L113+b333*L333));
  ka3(i)=k3(i)*km;
  ke3(i)=ka3(i)*fr3;
  
  ktim(i)=ke1(i)+ke2(i)+ke3(i);
  
end



%     plot(f,ke1)
% hold on;
%      plot(f,ke2)
%  hold on;
%    plot(f,ke3)
%  hold on;
  plot(f*100,ktim)
  hold on;
 
  
  
  
  
  fr21=0.70
  fr22=0.10
  fr23=0.20

f=0:0.001:0.1
% cs=1/3; %%%%
for i=1:length(f)
%   k11(i)=(2+ f(i)*(b11*(1-L11)*(1+cs)+ b33*(1-L33)*(1-cs) ))/(2-f(i)*(b11*L11*(1+cs)+b33*L33*(1-cs)));
 k1(i) =(3+f(i)*(2*b111*(1-L111)+b331*(1-L331)))/(3-f(i)*(2*b111*L111+b331*L331));
  ka1(i)=k1(i)*km;
  ke1(i)=ka1(i)*fr21;
  k2(i) =(3+f(i)*(2*b112*(1-L112)+b332*(1-L332)))/(3-f(i)*(2*b112*L112+b332*L332));
  ka2(i)=k2(i)*km;
  ke2(i)=ka2(i)*fr22;
   k3(i) =(3+f(i)*(2*b113*(1-L113)+b333*(1-L333)))/(3-f(i)*(2*b113*L113+b333*L333));
  ka3(i)=k3(i)*km;
  ke3(i)=ka3(i)*fr23;
  
  ktim2(i)=ke1(i)+ke2(i)+ke3(i);
  
end

    plot(f*100,ktim2)
    hold on;

     fr31=0.65
     fr32=0.10
     fr33=0.25

f=0:0.001:0.1
% cs=1/3; %%%%
for i=1:length(f)
%   k11(i)=(2+ f(i)*(b11*(1-L11)*(1+cs)+ b33*(1-L33)*(1-cs) ))/(2-f(i)*(b11*L11*(1+cs)+b33*L33*(1-cs)));
 k1(i) =(3+f(i)*(2*b111*(1-L111)+b331*(1-L331)))/(3-f(i)*(2*b111*L111+b331*L331));
  ka1(i)=k1(i)*km;
  ke1(i)=ka1(i)*fr31;
  k2(i) =(3+f(i)*(2*b112*(1-L112)+b332*(1-L332)))/(3-f(i)*(2*b112*L112+b332*L332));
  ka2(i)=k2(i)*km;
  ke2(i)=ka2(i)*fr32;
   k3(i) =(3+f(i)*(2*b113*(1-L113)+b333*(1-L333)))/(3-f(i)*(2*b113*L113+b333*L333));
  ka3(i)=k3(i)*km;
  ke3(i)=ka3(i)*fr33;
  
  ktim3(i)=ke1(i)+ke2(i)+ke3(i);
  
end

    plot(f*100,ktim3)

 xlabel(' Volume fraction f'), ylabel('Thermal conductivity of composite kc (W/mK)'),
% title('thermal conductivity vs. f=volume fraction for graphene-MLG-epoxy compodite: differnt sample ')
f1=f*100;
f2=f1';
iktim=ktim';
iktim2=ktim2';
iktim3=ktim3';
