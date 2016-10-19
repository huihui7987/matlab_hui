clc;
clear;


FWHM=62e-12;            %高斯信号FWHM宽度，为50ps
time_window=5000*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=8001;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

% n1=1.5;
% n2=0.66;
% n3=1.5;
%figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
%plot(t*1e+9,gauss_time,'linewidth',2.5);
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
%plot(double_f*1e-9,(gauss_spec),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('gauss_spec');hold

kk1 = 1./((1i*2*pi*double_f*1e-9)+0.016);
kk2 = 1./((1i*2*pi*double_f*1e-9)+0.025);
kk3 = 1./((1i*2*pi*double_f*1e-9)+0.035);
kk4 = 1./((1i*2*pi*double_f*1e-9)+0.045);
kk5 = 1./((1i*2*pi*double_f*1e-9)+0.055);
kk6 = 1./((1i*2*pi*double_f*1e-9)+0.065);
% figure;
% plot(double_f*1e-9,kk1,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('kk1-Intensity Transmission H1');hold on;
% figure;
hh1 = gauss_spec .* kk1;
% plot(double_f*1e-9,abs(hh),'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
figure;
hht1 = ifft(hh1);
hh2 = gauss_spec .* kk2;
hht2 = ifft(hh2);
hh3 = gauss_spec .* kk3;
hht3 = ifft(hh3);
hh4 = gauss_spec .* kk4;
hht4 = ifft(hh4);
hh5 = gauss_spec .* kk5;
hht5 = ifft(hh5);
hh6 = gauss_spec .* kk6;
hht6 = ifft(hh6);
plot(t*1e+9,abs(hht1)/(8.306e-6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht2)/(8.306e-6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht3)/(8.306e-6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht4)/(8.306e-6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht5)/(8.306e-6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht6)/(8.306e-6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
legend('k=0.016','k=0.025','k=0.035','k=0.045','k=0.055','k=0.065')



