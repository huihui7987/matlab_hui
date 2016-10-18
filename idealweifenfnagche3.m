clc;
clear;


FWHM=800e-12;            %高斯信号FWHM宽度，为50ps
time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=401;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

% n1=1.5;
% n2=0.66;
% n3=1.5;
figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
plot(t*1e+9,gauss_time,'linewidth',2.5);
% %以上，画图不美观，原因取点数太少
% xlabel('Time/ns');
% ylabel('Amplitude/V');
% title('Gauss pulse');

%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行移位操作。
%gauss_spec=fftshift(gauss_time); 
gauss_spec=gauss_spec/Ns;   %求实际的幅度值；归一化？
df=1/time_window;               %频率分辨率
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %双边频谱对应的频点
plot(double_f*1e-9,(gauss_spec),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('gauss_spec');hold

kk1 = 1./((1i*2*pi*double_f*1e-9)+0.048);
kk2 = 1./((1i*2*pi*double_f*1e-9)+0.028);
figure;
plot(double_f*1e-9,kk1,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('kk1-Intensity Transmission H1');hold on;
figure;
hh = gauss_spec .* kk1;
plot(double_f*1e-9,abs(hh),'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
figure;
hht = ifft(hh);
plot(t*1e+9,abs(hht)/0.001953,'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');




