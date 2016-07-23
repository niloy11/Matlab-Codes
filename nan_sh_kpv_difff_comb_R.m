clc;
close all;
clear all;
sigbA = 3.7*10^9 ;
Rb=[ 20 ]*10^-9;
%1 5 10 20
for l=1:length(Rb)
%%%%%%%%%%%
A1= 500 *10^-9;
L1= 300*10^-9;%%%%%%%%%%
% Lef1=(A* L1)/(A+ L1);
% v=2*10^4;
% C= 700*2200
% kp1 = 1/3 *Lef1* v* C
 kp1 =(sigbA )/(1/(L1)+2/(pi *A1));

%kp= 2000;
km=0.1;
H1=3*.35*10^-9;%%%%%%%%%%%
% H=.35*10^-9
a11=L1;
a21=L1;
a31=H1;

% kp1=550%%%%%%%%%%
%km=.1;



p1=H1/L1

% ak1=Rb*km
alp1(l)=Rb(l)*km/H1


L111=(p1^2/(2*(p1^2-1)))+(p1/(2*(1-p1^2)^(3/2)))*acos(p1);
gama1(l)=(1+2*p1)*alp1(l);

L221=L111;
L331=1-2*L111;

kc111(l)=kp1/(1+gama1(l)*L111*kp1/km);
kc221(l)=kp1/(1+gama1(l)*L221*kp1/km);
kc331(l)=kp1/(1+gama1(l)*L331*kp1/km);

% kc111=kp1 +gama1*(1-L111)*km;
% kc221=kp1 +gama1*(1-L221)*km;
% kc331=kp1 +gama1*(1-L331)*km;

b111(l)=(kc111(l)-km)/(km +L111*(kc111(l)-km));
b221(l)=(kc221(l)-km)/(km +L221*(kc221(l)-km));
b331(l)=(kc331(l)-km)/(km +L331*(kc331(l)-km));


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


% ak2=Rb*km
alp2(l)=Rb(l)*km/H2;

L112=(p2^2/(2*(p2^2-1)))+(p2/(2*(1-p2^2)^(3/2)))*acos(p2);
gama2(l)=(1+2*p2)*alp2(l);

L222=L112;
L332=1-2*L112;

kc112(l)=kp2/(1+gama2(l)*L112*kp2/km);
kc222(l)=kp2/(1+gama2(l)*L222*kp2/km);
kc332(l)=kp2/(1+gama2(l)*L332*kp2/km);

% kc112=kp2 +gama2*(1-L112)*km;
% kc222=kp2 +gama2*(1-L222)*km;
% kc332=kp2 +gama2*(1-L332)*km;

b112(l)=(kc112(l)-km)/(km +L112*(kc112(l)-km));
b222(l)=(kc222(l)-km)/(km +L222*(kc222(l)-km));
b332(l)=(kc332(l)-km)/(km +L332*(kc332(l)-km));




%%%%%%%%%%%
L3= 5000*10^-9;%%%%%%%%%%
H3=3*0.35*10^-9;
% Lef3=(A* L3)/(A+ L3);

   kp3=2200

a13=L3;
a23=L3;
a33=H3;

% kp3=2000%%%%%%%%%%%%%


p3=H3/L3

% ak3=Rb*km
alp3(l)=Rb(l)*km/H3


L113=(p3^2/(2*(p3^2-1)))+(p3/(2*(1-p3^2)^(3/2)))*acos(p3)
gama3(l)=(1+2*p3)*alp3(l);

L223=L113;
L333=1-2*L113;

kc113(l)=kp3/(1+gama3(l)*L113*kp3/km);
kc223(l)=kp3/(1+gama3(l)*L223*kp3/km);
kc333(l)=kp3/(1+gama3(l)*L333*kp3/km);

% kc113=kp3 +gama3*(1-L113)*km;
% kc223=kp3 +gama3*(1-L223)*km;
% kc333=kp3 +gama3*(1-L333)*km;

b113(l)=(kc113(l)-km)/(km +L113*(kc113(l)-km));
b223(l)=(kc223(l)-km)/(km +L223*(kc223(l)-km));
b333(l)=(kc333(l)-km)/(km +L333*(kc333(l)-km));



fr1=0.7
fr2=0.1
fr3=0.20





f=0:0.001:0.1
% cs=1/3; %%%%
for i=1:length(f)
%   k11(i)=(2+ f(i)*(b11*(1-L11)*(1+cs)+ b33*(1-L33)*(1-cs) ))/(2-f(i)*(b11*L11*(1+cs)+b33*L33*(1-cs)));
 k1(i) =(3+f(i)*(2*b111(l)*(1-L111)+b331(l)*(1-L331)))/(3-f(i)*(2*b111(l)*L111+b331(l)*L331));
  ka1(i)=k1(i)*km;
  ke1(i)=ka1(i)*fr1;
  k2(i) =(3+f(i)*(2*b112(l)*(1-L112)+b332(l)*(1-L332)))/(3-f(i)*(2*b112(l)*L112+b332(l)*L332));
  ka2(i)=k2(i)*km;
  ke2(i)=ka2(i)*fr2;
   k3(i) =(3+f(i)*(2*b113(l)*(1-L113)+b333(l)*(1-L333)))/(3-f(i)*(2*b113(l)*L113+b333(l)*L333));
  ka3(i)=k3(i)*km;
  ke3(i)=ka3(i)*fr3;
  
  ktim(i)=ke1(i)+ke2(i)+ke3(i);
  
end



%     plot(f,ke1)
% hold on;
% %     plot(f,ke2)
%  hold on;
% %   plot(f,ke3)
%  hold on;
  plot(f*100,ktim)
  hold on;
  xlabel(' volume fraction f'), ylabel('thermal conductivity of composite kc (W/mK)'),
% title('thermal conductivity vs. f=volume fraction for composite of graphene-MLG with various base TIM  ')

end
f2=(f*100)';
ktim1=ktim';
  