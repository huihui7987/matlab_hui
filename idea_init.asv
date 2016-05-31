clear all
R=50e-6;
lamda=1439.2e-9:1e-12:1439.8e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
neff=3.17995709;
%neff=3.17996;
r=0.83;

FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=601;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
%dt=1e-12;
t=0:dt:time_window;     %采样时间

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');
%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行移位操作。
gauss_spec=gauss_spec/Ns;   %求实际的幅度值；

df=1/time_window;               %频率分辨率
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %双边频谱对应的频点,频域横坐标



L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
taoa=0.83;
Ta1= exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));
%Ta1= (taoa-r)./(1-r.*taoa);
T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r<=taoa
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end
%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))
taob=0.801;
Ta2= exp(1i*(pi+phi)).*(taob-r.*exp(-1i*phi))./(1-r.*taob.*exp(1i*phi));
T2= (abs(Ta2)).^2;%~~~~-
PHI2 = angle(Ta2);%~~

taoc=0.846;
Ta3= exp(1i*(pi+phi)).*(taoc-r.*exp(-1i*phi))./(1-r.*taoc.*exp(1i*phi));
T3= (abs(Ta3)).^2;%~~~~-
PHI3 = angle(Ta3);%~~
if r<=taoa
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end


ffff =Ta2 .* gauss_spec;


subplot(1,2,1);
plot(v,T1,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
legend('Critical-coupled','Under-coupled','Over-coupled')
subplot(1,2,2);
plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
legend('Critical-coupled','Under-coupled','Over-coupled')
%subplot(1,3,3);plot(v,y3,'linewidth',2);title('n=1.5');
figure;
%subplot(1,3,3);
plot(double_f,abs(ffff).^2/3.5/1e-6/0.9397,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('res');hold on;

figure;
subplot(1,3,1);
plot(v,Ta2,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Ta1');hold on;
subplot(1,3,2);
plot(v,gauss_spec,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('gauss_spec');hold on;


