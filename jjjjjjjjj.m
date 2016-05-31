%==========================================================================
%Name:      spectrum_analysis.m
%Desc:      �Ը�˹�ź�Ϊ���������Ƶ�ס�˫�߹����ס����߹����ס�˫�߹������ܶȡ�
%           ���߹������ܶȣ������˹�źŵİ벨ȫ��FWHM=50ps�����ĵ�λ��2.5ns����
%=========================================================================
clc;
clear;
FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=50*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=1024;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-1.25e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');
%===========���¼���˫���ס�˫�߹����ס�˫�߹������ܶ�=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %����Ҷ�任�����ҽ���fftshift��λ������
gauss_spec=gauss_spec/Ns;   %��ʵ�ʵķ���ֵ��
df=1/time_window;               %Ƶ�ʷֱ���
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %˫��Ƶ�׶�Ӧ��Ƶ��


figure; %������
plot(double_f*1e-9,abs(gauss_spec),'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Amplitude/V');
title('double Amplitude spectrum');


figure; %��λ��
plot(double_f*1e-9,angle(gauss_spec),'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Phase/rad');
title('double Phase spectrum');


figure; %������
double_power_spec_W=abs(gauss_spec).^2;                 %˫�߹����ף���λW��
double_power_spec_mW=double_power_spec_W*1e+3;          %˫�߹����ף���λmW��
double_power_spec_dBm=10*log10(double_power_spec_mW);   %˫�߹����ף���λdBm��
plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');


figure; %�������ܶ�
double_power_specD_W=abs(gauss_spec).^2/(df);       %˫�߹������ܶ�,��λW/Hz
double_power_specD_mW=double_power_specD_W*1e+3;    %˫�߹������ܶ�,��λmW/Hz
double_power_specD_dBm=10*log10(double_power_specD_mW);%˫�߹������ܶ�,��λdBm/Hz
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