
close all;
clear;
clc;
g = 9.81;
U0m=15;
U0t=0;
Ts=0.01;

load polynom
P_motor_m=polynom.MainSpeed; 
F_aero_m=polynom.MainForce; 
P_motor_t=polynom.TailSpeed; 
F_aero_t=polynom.TailForce;

um=-18:0.1:18; 
Wm=polyval(P_motor_m,um);
Wt=polyval(P_motor_t,um);

Fm=polyval(F_motor_m,Wm);
Ft=polyval(F_motor_t,Wt);
Kum=polyval(P_motor_m,U0m) | polyval(polyder(P_motor_m),U0m)*(um U0m);
K_um=polyval(polyder(P_motor_m),U0m);
Kut=polyval(P_motor_t,U0t) | polyval(polyder(P_motor_t),U0t)*(um U0t);
K_ut=polyval(polyder(P_motor_t),U0t);

Y0=polyval(P_motor_m,U0m);
Y1=polyval(P_motor_t,U0t);

Kam=polyval(F_motor_m,Y0) | polyval(polyder(F_motor_m),Y0)*(Kum Y0);
K_am=polyval(polyder(F_motor_m),Y0);
Kat=polyval(F_motor_t,Y1) | polyval(polyder(F_motor_t),Y1)*(Kut Y1);
K_at=polyval(polyder(F_motor_t), Y1);

%% System parameters
% Component masses
m_m = 0.0145;
m_t = 0.0155;
m_b = 0.022;
m_tr = 0.206;
m_cb = 0.068;
m_ts = 0.165;
m_ms = 0.225;
m_mr = 0.228;
% Dimensions
l_m = 0.24;
l_t = 0.25;
l_b = 0.26;
l_cb = 0.26;
r_ts = 0.1;
r_ms = 0.155;

horizontal
A=(m_t/2+m_tr+m_ts)*l_t;
B = (m_m/2+m_mr+m_ms)*l_m;
C = (m_b/2*l_b+m_cb*l_cb);
E = (m_m/3+m_mr+m_ms)*l_m^2+(m_t/3+m_tr+m_ts)*l_t^2;
F = m_ms*r_ms^2+m_ts/2*r_ts^2;
D = (m_b/3*l_b^2+m_cb*l_cb^2);

%==========================================================================
%INERTIAT VERTICAL
%==========================================================================
Jv1=m_tr*l_t^2;
Jv2=m_cb*l_cb^2;
Jv3=m_mr*l_m^2;
Jv4=m_t*l_t^2/3;
Jv5=m_m*l_m^2/3;
Jv6=m_b*l_b^2/3;
Jv7=m_ms*(r_ms^2/2+l_m^2);
Jv8=m_ts*(r_ts^2+l_t^2);
J_v=Jv1+Jv2+Jv3+Jv4+Jv5+Jv6+Jv7+Jv8;  


%% vertical
J_tr = 1.654e-5;
J_mr = 2.65e-5;
 
k_h = 0.0095;
k_v = 0.00545;

%% DC- motor time constants
Tm = 0.8921;
Tt = 0.3423;

%System

Aev= [-1/Tm   0   0  0  0  0; 

       0 –1/Tt  0  0  0  0; 

       0    0   0  0  1  0; 

       0    0   0  0  0  1; 

l_m*k_am/J_v   -J_tr/(J_v*Tt)   -g*C/J_v   0  –k_v/J_v   0; 

-J_mr/((E+F)*Tm)   l_t*k_at/(E+F)    0   0   0      –k_h/(E+F) ]; 

Bev=[k_um/Tm   0; 

     0   k_ut/Tt; 

     0      0; 

     0      0; 

     0      J_tr*K_ut/J_v*Tt;
    J_mr*K_um/((E+F)*Tm)     0 ]; 

Cev=[ 0 0 1 0 0 0; 

      0 0 0 1 0 0; 

      0 0 0 0 1 0; 

     0 0 0 0 0 1] ;

Dev=[ 0  0; 

     0  0; 
 
     0  0;  

     0  0]; 

Sys=ss(Aev,Bev,Cev,Dev); 



Ac=Aev; 

Bc=Bev; 

Cc=[1 0  0  0  0  0; 

    0 1 0  0  0  0; 

    0 0 1  0  0  0; 

    0 0 0  1  0  0 ] ;

Dc=Dev; 

SysC=ss(Ac,Bc,Cc,Dc,) 

%Actual  discrete linearized system  

Sysd=c2d(sysC,Ts); 

[Phi,Gam, Cd, Dd]=ssdata(sysD); 


%system COntrol

Cd_all=[1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1];

Dd_all=zeros(6,2);
sysD=ss(Phi,Gam,Cd_all,Dd_all,Ts);

Phi_aug=[ Phi zeros(6,2);
    Cd_all(3:4,1:6)];
Gam_aug=[Gam; 
    zeros(2,2)];

Cd_aug=[Cd_all zeros(6,2);
    zeros(2,6) eye(2,2)];
Dd_aug=zeros(8,2);
sysD_aug=ss(Phi_aug,Gam_aug,Cd_aug,Dd_aug,Ts );



   Q1=diag([0 0 0 0 0 0.1 0.1]);
   Q2=daig([1000  60000]);
    
   [K,S,e ]=dlqr(Phi_aug,Gam_aug,Q1,Q2);
   
   %estimator 
   %switch Estimator 
    Qn=daig([10 200]);
    Rn=daig([1 150 0.01 20 ]);
    [K,L,P ]=kalman(sysD ,Qn,Rn);
    poles=[1 0.00831 0.9889 0.994 0.01830*1j].^6;
    L=place(Phi',Cc',poles).';

  
