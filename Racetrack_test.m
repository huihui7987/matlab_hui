clear all
R=20e-6;
lamda=1562.7e-9:1e-12:1563.1e-9;
v=(3e8./lamda)-(3e8./1562.865e-9);
%neff=3.179992;
neff=2.5;
r=0.98;
%yt=0.999;
Lc = R;
L = 2*pi*R+2*Lc;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi);
taoa=0.98;

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

taob=0.97;
Ta2= (r-taob.*p)./(1-r.*taob.*p);
T2= (abs(Ta2)).^2;
PHI2 = angle(Ta2);
if r<=taob
    PHI2 = PHI2+(PHI2<0)*2*pi ;
end

taoc=0.99;
Ta3= (r-taoc.*p)./(1-r.*taoc.*p);
T3= (abs(Ta3)).^2;
PHI3 = angle(Ta3);
if r<=taoc
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end

%subplot(1,2,1);
plot(v,T1,'r','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Intensity Transmission');hold on;
plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Intensity Transmission');hold on;
plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Intensity Transmission');hold on;
legend('Critical-coupled','Under-coupled','Over-coupled')
figure;
%subplot(1,2,2);
plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Phase Response');hold on;
plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Phase Response');hold on;
plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Phase Response');hold on;
legend('Critical-coupled','Under-coupled','Over-coupled')
figure;

plot(v,TT1,'g','linewidth',2); xlabel('Frequency(Hz£©');ylabel('Tran_Response');hold on;
