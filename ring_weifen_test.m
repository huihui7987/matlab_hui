clc;
clear;


FWHM=62e-12;            %高斯信号FWHM宽度，为50ps
time_window=3000*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=1001;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

% n1=1.5;
% n2=0.66;
% n3=1.5;
%figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
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

%plot(double_f*1e-9,gauss_spec,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('gau');hold on;
%%

R=50e-6;
lamda=1439e-9:1e-12:1440e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
%neff=3.179962;
%neff=3.289700581;
%neff = 3.28489059783;
neff = 3.224659641;

%yt=0.999;

L = 2*pi*R;
%phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
phi = mod(L*neff./lamda*pi*2,6*pi);%~~~~-
p=exp(1i*phi/2);
taoa=0.9999;

%figure;
%Tdaa= taoa.*(r-1)./(1-r.*taoa.*p) ;
r1 = 0.895;
r2 = 0.999;
k1 = sqrt(1-r1^2);
k2 = sqrt(1-r2^2);
det = 5e-9;
Tdaa = (-taoa*k1*k2*p)./(1-taoa*r1*r2*p.^2);

H_ring_res1 = gauss_spec .* Tdaa;

%plot(double_f*1e-9,abs(Tdaa),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('传输函数');hold on;

ring_gauss_diff_power_spec=(abs(H_ring_res1)).^2;
%figure;
%plot(double_f*1e-9,abs(H_ring_res1),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('频域结果');hold on;
figure;
hht = abs(ifft((H_ring_res1)))/9.566e-8/2.878/0.3475/0.765;
plot(t*1e+9,hht,'linewidth',2.5);xlabel('Time(ps）');ylabel('Intensity(a.u.)');
figure;
ff = angle(Tdaa);
plot(double_f*1e-9,ff,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('Tdaa');title('gau');hold on;

