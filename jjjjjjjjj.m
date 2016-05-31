%==========================================================================
%Name:      spectrum_analysis.m
%Desc:      以高斯信号为例，求解其频谱、双边功率谱、单边功率谱、双边功率谱密度、
%           单边功率谱密度，这里高斯信号的半波全宽FWHM=50ps，中心点位于2.5ns处。
%=========================================================================
clc;
clear;
FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
time_window=50*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=1024;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-1.25e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');
%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行fftshift移位操作。
gauss_spec=gauss_spec/Ns;   %求实际的幅度值；
df=1/time_window;               %频率分辨率
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %双边频谱对应的频点


figure; %幅度谱
plot(double_f*1e-9,abs(gauss_spec),'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Amplitude/V');
title('double Amplitude spectrum');


figure; %相位谱
plot(double_f*1e-9,angle(gauss_spec),'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Phase/rad');
title('double Phase spectrum');


figure; %功率谱
double_power_spec_W=abs(gauss_spec).^2;                 %双边功率谱，单位W；
double_power_spec_mW=double_power_spec_W*1e+3;          %双边功率谱，单位mW；
double_power_spec_dBm=10*log10(double_power_spec_mW);   %双边功率谱，单位dBm；
plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');


figure; %功率谱密度
double_power_specD_W=abs(gauss_spec).^2/(df);       %双边功率谱密度,单位W/Hz
double_power_specD_mW=double_power_specD_W*1e+3;    %双边功率谱密度,单位mW/Hz
double_power_specD_dBm=10*log10(double_power_specD_mW);%双边功率谱密度,单位dBm/Hz
plot(double_f*1e-9,double_power_specD_dBm,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/(dBm/Hz)');
title('double power spectrum Density');


%%
racetrack
%%
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