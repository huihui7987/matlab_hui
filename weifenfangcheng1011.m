clear all
R=50e-6;
lamda=1439.2e-9:1e-12:1439.8e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
neff=3.179962;
%neff=2.5;
r=0.83;
%yt=0.999;
Lc = R;
L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi);
taoa=0.83;

%Ta1= exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));
Ta1= (r-taoa.*p)./(1-r.*taoa.*p);
%Ta1= exp(1i*(pi+phi)).*(taoa*(yt)-r.*exp(-1i*phi))./(1-r.*taoa*(yt).*exp(1i*phi));

%Ta1=r*((1-yt*p)/(1-r^2*yt*p))

TT1=abs(Ta1);
T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r<=taoa
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end
%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))

taob=0.801;
Ta2= (r-taob.*p)./(1-r.*taob.*p);
T2= (abs(Ta2)).^2;
PHI2 = angle(Ta2);
if r<=taob
    PHI2 = PHI2+(PHI2<0)*2*pi ;
end

taoc=0.846;
Ta3= (r-taoc.*p)./(1-r.*taoc.*p);
T3= (abs(Ta3)).^2;
PHI3 = angle(Ta3);
if r<=taoc
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end

%subplot(1,2,1);
% plot(v,T1,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% legend('Critical-coupled','Under-coupled','Over-coupled')
% figure;
% %subplot(1,2,2);
% plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
% plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
% plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
% legend('Critical-coupled','Under-coupled','Over-coupled')

%%����һ�����ߣ�Ϊ΢���Ĵ��亯�����ı书�ʿ��Ըı䴫�亯��
figure;
Tda= taoa.*(r-1).^2./(1-r.*taoa.*p).*2 ;
Tdb= taob.*(r-1).^2./(1-r.*taob.*p).*2 ;
Tdc= taoc.*(r-1).^2./(1-r.*taoc.*p).*2 ;
%ttt = (abs((1-r).^2.*taob.^0.5.*p./(1-r.*taob.*p))).^2;
plot(v,Tda/0.1642,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
plot(v,Tdb/0.1642,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
plot(v,Tdc/0.1642,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
legend('pump=xx','pump=x','pump=xxx')
%plot(v,TTT,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Tran_Response');hold on;
%%����͸����
% n1 = 0.3;
% n2 = 0.6;
% n3 = 1;
% %H1_iii = abs(1./(1i*2*pi*v+0.5)).^0.6;
% H1_ideal=1./(abs(1i*2*pi*v*1e-9)).^(2*n1);
% H2_ideal=(1./(abs(1i*2*pi*v*1e-9)).^(2*n2));
% H3_ideal=(1./(abs(1i*2*pi*v*1e-9)).^(2*n3));
% figure;
% subplot(1,3,1);
% plot(v,H1_ideal,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
% subplot(1,3,2);
% plot(v,H2_ideal,'g','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H2');hold on;
% subplot(1,3,3);
% plot(v,H3_ideal,'b','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H3');hold on;
figure;
kk1 = 1./(abs(2*pi*v*1e-9)+0.33).^0.2;
kk2 = 1./(abs(2*pi*v*1e-9)+0.48).^0.2;
kk3 = 1./(abs(2*pi*v*1e-9)+0.57).^0.2;
plot(v,kk1/1.248,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
plot(v,kk2/1.248,'g','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
plot(v,kk3/1.248,'b','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
legend('pump=0.33','pump=.48','pump=.57')
