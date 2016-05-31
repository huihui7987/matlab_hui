clc;
clear;
FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=2048;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

n1=1.5;
n2=1;
n3=1.5;

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');

%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行移位操作。
gauss_spec=gauss_spec/Ns;   %求实际的幅度值；归一化？
df=1/time_window;               %频率分辨率
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %双边频谱对应的频点

% figure;
% 
% plot(double_f*1e-9,gauss_spec,'linewidth',2.5);
% figure; %幅度谱
% plot(double_f*1e-9,abs(gauss_spec),'linewidth',2.5);
% xlabel('Frequency/GHz');
% ylabel('Amplitude/V');
% title('double Amplitude spectrum');
% 
% 
% figure; %相位谱
% plot(double_f*1e-9,angle(gauss_spec),'linewidth',2.5);
% xlabel('Frequency/GHz');
% ylabel('Phase/rad');
% title('double Phase spectrum');

figure; %高斯脉冲功率谱
double_power_spec_W=abs(gauss_spec).^2;                 %双边功率谱，单位W；
double_power_spec_mW=double_power_spec_W*1e+3;          %双边功率谱，单位mW；
double_power_spec_dBm=10*log10(double_power_spec_mW);   %双边功率谱，单位dBm；
plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');
figure;


%%%%%%%%%%%%%%%%%
%%理想微分器传输函数
%%%实际证明该函数可用，为复数形式，matlab画图只能画出其实部，幅度与相位可分别利用abs()与angle()查看；0512
%%%%%%%%%%%%%%%%%
%传输函数
H1=(1i*2*pi*double_f*1e-9).^(n1);
%%%分段，暂时不知道错误在哪里，目前确定不对；0512
%H1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
idea_Intensity_Transmission = mapminmax((abs(1i*2*pi*double_f*1e-9)).^(2*n1),0,1);
% idea_power_spec = (abs(H1)).^2;
 idea_phase_spec = angle(H1);
 
%%%%%作图%%%%%%
%%%%%%%%%%%%%%%
 %理想传输函数
 subplot(1,3,1)
 plot(double_f*1e-9,abs(H1),'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('idea_tran_spec');title('理想传输函数');hold on;
 %理想传输函数相位
 subplot(1,3,2)
 plot(double_f*1e-9,idea_phase_spec,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('idea_phase_spec');title('理想相位变化');hold on;
 %理想传输函数透射谱
subplot(1,3,3)
plot(double_f*1e-9,idea_Intensity_Transmission,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');title('理想透射谱');hold on;


%%%%%高斯函数理想微分之后的功率谱%%%%%
figure;
%理想微分结果
idea_gauss_diff=gauss_spec.*H1;
%理想微分功率谱
idea_gauss_diff_power_spec=(abs(idea_gauss_diff)).^2;
subplot(1,2,1)
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');title('理想微分功率谱');hold on;
subplot(1,2,2)
%%%微分输出结果
plot(double_f*1e-9,abs(idea_gauss_diff),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('idea_gauss_diff(a.u.)');hold on;


%%微环微分器传输函数
%n=1.5
figure;
t0=-3*1e-4;
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
H3=(-1i*2*pi*double_f*1e-9).^(n1);
Hr1= -exp(-1i*double_f*1e-9*t0).*H3;
Hr1_tran = abs(Hr1).^2;
Hr1_Pi = angle(Hr1);

%微环传输函数
subplot(1,3,1)
plot(double_f*1e-9,Hr1,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('MRR_tran_spec');title('微环传输函数');hold on;
subplot(1,3,2)
plot(double_f*1e-9,Hr1_Pi,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('MRR_tran_spec');hold on;
subplot(1,3,3)
plot(double_f*1e-9,Hr1_tran,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('MRR_tran_spec');hold on;


%%%%%%%高斯函数经过微环微分器之后%%%%%%%%%%
figure;
%t0=-3*1e-4;
%H1=(1i*2*pi*double_f*1e-9).^(n1);
%H3=(-1i*2*pi*double_f*1e-9).^(n1);
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
%传输函数
%Hr1= -exp(-1i*double_f*1e-9*t0).*H3;
%微分结果
MRR_gauss_diff = gauss_spec.*Hr1;
%微环微分之后功率谱
MRR_gauss_diff_power_spec=abs(MRR_gauss_diff).^2;

subplot(1,2,1)
plot(double_f*1e-9,MRR_gauss_diff_power_spec,'g','linewidth',2.5);
subplot(1,2,2)
plot(double_f*1e-9,MRR_gauss_diff,'g','linewidth',2.5);


%%%%%%%%%%test%%%%%%%%%%
% 
% R=20e-6;
% %neff=3.179992;
% neff=2.5;
% r=0.98;
% %yt=0.999;
% Lc = R;
% L = 2*pi*R+2*Lc;
% phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
% p=exp(1i*phi);
% taoa=0.98;
% 
% %Ta1= exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));
% Ta1= (r-taoa.*p)./(1-r.*taoa.*p);
% %Ta1= exp(1i*(pi+phi)).*(taoa*(yt)-r.*exp(-1i*phi))./(1-r.*taoa*(yt).*exp(1i*phi));
% TT1=abs(Ta1);
% %Ta1=r*((1-yt*p)/(1-r^2*yt*p))
% figure;
% 
% plot(v,TT1,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Tran_Response');hold on;






