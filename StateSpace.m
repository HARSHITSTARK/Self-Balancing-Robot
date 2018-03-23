close all
clear all
clc

%{
A = [0 1 0 0; 62.0193 -44.5897 0 -2123.32; 0 0 0 1; 6.09908 -10.1911 0 -485.289]; 
B = [0; -90.0275; 0; -20.5759]; 
C = [1 0 0 0; 0 0 1 0];          
D = zeros(size(C,1),size(B,2)); 
%}


kt=.3233;
kb=.4953;
R=5.2628;
L=.11;
mp=.3438;
mw=.0352;
rw=.021;
ip=.0042;
icmw=7.832e-6;
g=9.8;

Mt=icmw*(ip+L.^2*mp)+(L.^2*mp*mw+ip*(mp+mw))*(rw).^2;

A=[0 1 0 0;
    (g*L*mp*(icmw+(mp+mw)*(rw).^2))/Mt (-kb*kt*(icmw+rw*(mw*rw+mp*(L+rw))))/(R*Mt) 0 -(kb*kt*(icmw+rw*(mw*rw+mp*(L+rw))))/(R*rw*Mt);
    0 0 0 1;
    (g*L.^2*(mp).^2*(rw).^2)/Mt (-kb*kt*rw*(ip+L*mp*(L+rw)))/(R*Mt) 0 -(kb*kt*(ip+L*mp*(L+rw)))/(R*Mt)];
 
B=[0; (-kt*(icmw+rw*(mw*rw+mp*(L+rw))))/(R*Mt); 0; (-kt*rw*(ip+L*mp*(L+rw)))/(R*Mt)];

C = [1 0 0 0;
     0 0 1 0];
 
D = [0;
     0];

syms s;
pole = eig(A); 
Sys1 = ss(A, B, C, D); 

Sys1_char = charpoly(A,s);

[Sys1num, Sys1den]=ss2tf(A,B,C,D);



Q = C'*C;
Q(1,1) = 30000;
Q(3,3) = 2000;
R = 1;
K = lqr(A,B,Q,R);
Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

sys_lqr = ss(Ac,Bc,Cc,Dc);

t = 0:0.01:5;
[y,t,x]=step(sys_lqr, t(end));

yyaxis left
plot(t,y(:,1));
ylabel('Angle From Vertical Position (radians)')
xlabel('time')

yyaxis right
plot(t,y(:,2));
ylabel('Miseg Translational Position (meters)')
title('Step Response of Minseg Robot with LQR Control')