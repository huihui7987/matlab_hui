
clear all;
clc
c       = 3e8;%����

%%%%%%%%%%%%%%%%%%%%%%%%%
%E1%
%E2%���
%E3%��2*pi*R��
%E4%��0*pi*R��
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













