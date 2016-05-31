clear all;
clc
c       = 3e8;

%%%%%%%%%%%%%%%%%%%%%%%%%
E1=1;%输入%归一化
%E2%输出
%E3%环2*pi*R处
%E4%环0*pi*R处
%lambda%波长
%neff%有效折射率
%beta=2*pi*neff/lambda;%传输常数

%%%%%%%%%%%%%%%%%%%%%%%%%
%R%半径
%L=2*pi*R;%周长

tao=0.95;%amplitude transmission factor
F=[4.4,12,30,47,56];%pi/(1-r*tao)
PIlinshi=repmat(pi,1,5);
r=(1-PIlinshi./F)/tao;
t=sqrt(1-r.^2);

%phi=beta*L;%round-trip phase shift
%phi=-3*pi/2:pi/2000:3*pi/2;%为作图，横坐标
phi=-0.3*pi:pi/2000:0.3*pi;
%T% intensity transmission factor
%PHI%"phase transmission factor"
%F=pi/(1-r*tao);
%Q=neff*F*L/lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%START
for cnd=1:5
    r1=r(cnd);
    t1=t(cnd);
    
%     
%   E4=r1*E3+1i*t1*E1;
%     E3=tao*exp(1i*phi)*E4;
%     
    
% E2=r1*E1+1i*t1*E3;
% E4=r1*E3+1i*t1*E1;
% E3=tao*exp(1i*phi)*E4;
% 
% T=(abs(E2/E1))^2;
% PHI=angle(E2/E1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ta = exp(1i*(pi+phi)).*(tao-r.*exp(-1i*phi))./(1-r.*tao.*exp(1i*phi));
T = (abs(Ta)).^2;%~~~~-

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(2,1,1)
plot(phi,Ta)
% subplot(2,1,2)
% plot(phi,PHI)



end;









