clear all
R=50e-6;
lamda=1439.2e-9:1e-12:1439.8e-9;
v=(3e8./lamda)-(3e8./1439.5e-9)
neff=3.17996;
%neff=2.5688;
%%%%%%%根tao的大小没有关系，而是与tao距离的关系
tao = 0.85 ;
r1 = 0.85;
r2 = 0.845;
r3 = 0.854;

L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-


Ta1= exp(1i*(pi+phi)).*(tao-r1.*exp(-1i*phi))./(1-r1.*tao.*exp(1i*phi));
T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r1<=tao
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end


%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))

Ta2= exp(1i*(pi+phi)).*(tao-r2.*exp(-1i*phi))./(1-r2.*tao.*exp(1i*phi));
T2= (abs(Ta2)).^2;%~~~~-
PHI2 = angle(Ta2);%~~
if r2<=tao
    PHI2 = PHI2+(PHI2<0)*2*pi ;
end

Ta3= exp(1i*(pi+phi)).*(tao-r3.*exp(-1i*phi))./(1-r3.*tao.*exp(1i*phi));
T3= (abs(Ta3)).^2;%~~~~-
PHI3 = angle(Ta3);%~~
if r3<=tao
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end


subplot(1,2,1);
plot(v,T1,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
legend('r1','r2','r3')
subplot(1,2,2);
plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
legend('r1','r2','r3')
%subplot(1,3,3);plot(v,y3,'linewidth',2);title('n=1.5');
