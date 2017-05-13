clc;
clear;

FWHM=50e-12;            %高斯信号FWHM宽度，为62ps
time_window=20000*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=801;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

% n1=1.5;
% n2=0.66;
% n3=1.5;
% figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2);%.*exp(-1i*1439.5e-9*t); %高斯脉冲，中心位于2.5ns处。
% plot(t*1e+9,gauss_time,'linewidth',2.5);
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
figure;
plot(double_f*1e-9,(gauss_spec),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('gauss_spec');hold
%% 理想传输函数
kk1 = 1./((1i*2*pi*double_f*1e-9)+0.033);
kk2 = 1./((1i*2*pi*double_f*1e-9)+0.035);
kk3 = 1./((1i*2*pi*double_f*1e-9)+0.038);
kk4 = 1./((1i*2*pi*double_f*1e-9)+0.055);
kk5 = 1./((1i*2*pi*double_f*1e-9)+0.068);
kk6 = 1./((1i*2*pi*double_f*1e-9)+0.081);

figure;%理想相位变化
plot(double_f*1e-9,angle(kk1),'r','linewidth',2.5); xlabel('Frequency( GHz）');ylabel('Intensity(a.u.)');hold on;
figure;%理想透射谱
res_rr1 = abs(kk1.^2);
res_rr2 = abs(kk2.^2);
res_rr3 = abs(kk3.^2);
res_rr4 = abs(kk4.^2);
res_rr5 = abs(kk5.^2);
plot(double_f*1e-9,res_rr1,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Intensity(a.u.)');hold on;
plot(double_f*1e-9,res_rr2,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Intensity(a.u.)');hold on;
plot(double_f*1e-9,res_rr3,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Intensity(a.u.)');hold on;
plot(double_f*1e-9,res_rr4,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Intensity(a.u.)');hold on;
plot(double_f*1e-9,res_rr5,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Intensity(a.u.)');hold on;
legend('k=0.031','k=0.043','k=0.052','k=0.064','k=0.071')

hh1 = gauss_spec .* kk1;
% plot(double_f*1e-9,abs(hh),'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
figure;
hht1 = abs(ifft(hh1))/max(abs(ifft(hh1)));
hh2 = gauss_spec .* kk2;
hht2 = abs(ifft(hh2))/max(abs(ifft(hh2)));
hh3 = gauss_spec .* kk3;
hht3 = abs(ifft(hh3))/max(abs(ifft(hh3)));
hh4 = gauss_spec .* kk4;
hht4 = abs(ifft(hh4))/max(abs(ifft(hh4)));
hh5 = gauss_spec .* kk5;
hht5 = abs(ifft(hh5))/max(abs(ifft(hh5)));
hh6 = gauss_spec .* kk6;
hht6 = abs(ifft(hh6))/max(abs(ifft(hh6)));
plot(t*1e+9,abs(hht1),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht2),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht3),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht4),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht5),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
legend('k=0.031','k=0.043','k=0.052','k=0.064','k=0.071','k=0.083')
% figure;
% hf1 = gauss_spec .* abs(kk1);
% hf2 = abs(ifft(hf1));
% hf3 = gauss_spec .* abs(kk5);
% hf4 = abs(ifft(hf3));
% plot(t*1e+9,abs(hf2),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
% plot(t*1e+9,abs(hf4),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;


%% 微环谐振腔实现
R=50e-6;
lamda=1439.1e-9:1e-12:1439.9e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
%neff=3.179962;
neff=3.40906;

L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi/2);
taoa=0.99;

figure;
%Tdaa= taoa.*(r-1)./(1-r.*taoa.*p) ;
r1 = 0.99;
r2 = 0.99;
k1 = sqrt(1-r1^2);
k2 = sqrt(1-r2^2);

taob = 0.95;
taoc = 0.90;
taod = 0.85;

Tdaa = ((-taoa*k1*k2*p)./(1-taoa*r1*r2*p.^2));
Tdaa2 = ((-taob*k2*k2*p)./(1-taob*r1*r2*p.^2));
Tdaa3 = ((-taoc*k1*k2*p)./(1-taoc*r1*r2*p.^2));
Tdaa4 = ((-taod*k1*k2*p)./(1-taod*r1*r2*p.^2));

H_ring_res1 = gauss_spec .* Tdaa;
%semilogy
plot(double_f*1e-9,20*log10(abs(Tdaa)),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('传输函数功率谱');hold on;
plot(double_f*1e-9,20*log10(abs(Tdaa2)),'g','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('传输函数功率谱');hold on;
plot(double_f*1e-9,20*log10(abs(Tdaa3)),'b','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('传输函数功率谱');hold on;
plot(double_f*1e-9,20*log10(abs(Tdaa4)),'y','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('传输函数功率谱');hold on;
% figure;%理想相位变化
% PHI1 = angle(Tdaa);
% 
% plot(double_f*1e-9,PHI1,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Tdaa-Intensity Transmission H1');hold on;
% ring_gauss_diff_power_spec=(abs(H_ring_res1)).^2;
figure;
plot(double_f*1e-9,abs(H_ring_res1),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('频域结果');hold on;
figure;
hht = abs(ifft((H_ring_res1)))/max(abs(ifft((H_ring_res1))));
plot(t*1e+9,hht,'linewidth',2.5);

%%

figure;
plot((t)*1e+9,abs(hht1),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t)*1e+9,abs(hht2),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t)*1e+9,abs(hht3),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t)*1e+9,abs(hht4),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t)*1e+9,abs(hht5),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t)*1e+9,abs(hht6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t-0.000000001)*1e+9,hht,'r','linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
legend('k=0.038','k=0.054','k=0.070','k=0.084','k=0.2','k=0.3','ring-based')


figure;
plot((t)*1e+9,abs(hht6),'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
plot((t-0.000000005)*1e+9,hht,'r','linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');hold on;
legend('Calculated','Ring-based')







