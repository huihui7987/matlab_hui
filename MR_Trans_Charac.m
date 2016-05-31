
clear all;
clc
c       = 3e8;%光速

%%%%%%%%%%%%%%%%%%%%%%%%%
%E1%
%E2%输出
%E3%环2*pi*R处
%E4%环0*pi*R处
lambda%
neff%
beta=2*pi*neff/lambda;%

%%%%%%%%%%%%%%%%%%%%%%%%%
R%
L=2*pi*R;%

r%
t=sqrt(1-r^2);%
tao%amplitude transmission factor
phi=beta*L;%round-trip phase shift
%T% intensity transmission factor
%PHI%"phase transmission factor"
F=pi/(1-r*tao);
Q=neff*F*L/lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%START
E2=r*E1+1i*t*E3;
E4=r*E3+1i*t*E1;
E3=tao*exp(1i*phi)*E4;

T=(abs(E2/E1))^2;
PHI=angle(E2/E1);













