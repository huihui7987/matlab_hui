clear all
lamda=1562.7e-9:1e-12:1563.1e-9;
v=(3e8./lamda)-(3e8./1562.9e-9);

n1=1;
% n2=1;
% n3=1.5;
% p1=(n1*pi/2).*(v>0)+(-n1*pi/2).*(v<0);
% p2=(n2*pi/2).*(v>0)+(-n2*pi/2).*(v<0);
% p3=(n3*pi/2).*(v>0)+(-n3*pi/2).*(v<0);
% 
% plot(v,p1,'r','linewidth',2.5);hold on;
% plot(v,p2,'g','linewidth',2.5); hold on;
% plot(v,p3,'b','linewidth',2.5);hold on;
% figure;
H1=((1i*2*pi*v)).^(n1);
% H2=(1i*2*pi*v).^(n2);
% H3=(1i*2*pi*v).^(n3);
H11=abs(H1);
PI = angle(H1)
subplot(1,2,1)
plot(v,PI,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
subplot(1,2,2)
plot(v,H1,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% plot(v,H2,'g'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% plot(v,H3,'b'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% legend('0.5','1','1.5')

figure;
H11=((abs(2*pi*v)).*(n1)).*exp(1i*n1*pi/2).*(v>0)+((abs(2*pi*v)).*(n1)).*exp(-1i*n1*pi/2).*(v<0);

plot(v,H11,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;



% 
% clc;
% clear;
% FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
% time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
% Ns=2048;                %采样点
% dt=time_window/(Ns-1);  %采样时间间隔
% t=0:dt:time_window;     %采样时间
% 
% n1=1.5;
% n2=1;
% n3=1.5;
% 
% gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
% %plot(t*1e+9,gauss_time,'linewidth',2.5);
% xlabel('Time/ns');
% ylabel('Amplitude/V');
% title('Gauss pulse');
% 
% %===========以下计算双边谱、双边功率谱、双边功率谱密度=================
% gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行移位操作。
% gauss_spec=gauss_spec/Ns;   %求实际的幅度值；
% df=1/time_window;               %频率分辨率
% k=floor(-(Ns-1)/2:(Ns-1)/2);    
% % k=0:Ns-1;
% double_f=k*df;   %双边频谱对应的频点
% 
% % figure;
% % 
% % plot(double_f*1e-9,gauss_spec,'linewidth',2.5);
% % figure; %幅度谱
% % plot(double_f*1e-9,abs(gauss_spec),'linewidth',2.5);
% % xlabel('Frequency/GHz');
% % ylabel('Amplitude/V');
% % title('double Amplitude spectrum');
% % 
% % 
% % figure; %相位谱
% % plot(double_f*1e-9,angle(gauss_spec),'linewidth',2.5);
% % xlabel('Frequency/GHz');
% % ylabel('Phase/rad');
% % title('double Phase spectrum');
% 
% %figure; %高斯脉冲功率谱
% double_power_spec_W=abs(gauss_spec).^2;                 %双边功率谱，单位W；
% double_power_spec_mW=double_power_spec_W*1e+3;          %双边功率谱，单位mW；
% double_power_spec_dBm=10*log10(double_power_spec_mW);   %双边功率谱，单位dBm；
% %plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
% xlabel('Frequency/GHz');
% ylabel('Power/dBm');
% title('double Power spectrum');
% 
% %%%%%%%%%%%%%%%%%
% %%理想微分器传输函数
% %%%%%%%%%%%%%%%%%
% id_H1=(1i*2*pi*double_f*1e-9).^(n1);%%%总
% id_Hr1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
% subplot(1,2,1);
% plot(double_f*1e-9,id_H1,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% subplot(1,2,2)
% plot(double_f*1e-9,id_Hr1,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% 
% %%%微环中n=1.5
% 
% t0=-3*1e-4;
% Mr_H1=(exp(-1i*double_f*1e-9*t0)).*(-1i*2*pi*double_f*1e-9).^(n1);%%%总
% %%分
% H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
% Mr_Hr1= exp(-1i*double_f*1e-9*t0).*H3;
% 
% figure;
% subplot(1,2,1);
% plot(double_f*1e-9,Mr_H1,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% subplot(1,2,2);
% plot(double_f*1e-9,Mr_Hr1,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;






    
