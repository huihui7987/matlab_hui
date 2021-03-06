clc;
clear;
R=50e-6;
lamda=1439.2e-9:1e-12:1439.8e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
neff=3.179962;
%neff=2.5;
r=0.9999;
%yt=0.999;
Lc = R;
L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi);
taoa=0.83;

FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=601;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

% n1=1.5;
% n2=0.66;
% n3=1.5;
figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
plot(t*1e+9,gauss_time,'linewidth',2.5);
%以上，画图不美观，原因取点数太少
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');

%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行移位操作。
%gauss_spec=fftshift(gauss_time); 
gauss_spec=gauss_spec/Ns;   %求实际的幅度值；归一化？
df=1/time_window;               %频率分辨率
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %双边频谱对应的频点
figure; %高斯脉冲功率谱
double_power_spec_W=abs(gauss_spec).^2;                 %双边功率谱，单位W；
double_power_spec_mW=double_power_spec_W*1e+3;          %双边功率谱，单位mW；
double_power_spec_dBm=10*log10(double_power_spec_mW);   %双边功率谱，单位dBm；
plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

taob=0.801;
Ta2= (r-taob.*p)./(1-r.*taob.*p);
T2= (abs(Ta2)).^2;
PHI2 = angle(Ta2);
if r<=taob
    PHI2 = PHI2+(PHI2<0)*2*pi ;
end

taoc=0.846;
Ta3= (r-taoc.*p)./(1-r.*taoc.*p);
T3= (abs(Ta3)).^2;
PHI3 = angle(Ta3);
if r<=taoc
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end

%subplot(1,2,1);
% plot(v,T1,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
% legend('Critical-coupled','Under-coupled','Over-coupled')
% figure;
% %subplot(1,2,2);
% plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
% plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
% plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
% legend('Critical-coupled','Under-coupled','Over-coupled')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 画出一组曲线，为微环的传输函数，改变功率可以改变传输函数
figure;%透射谱，是传输函数的平方
Tda= taoa.*(r-1).^2./(1-r.*taoa.*p).*2 ;
Tdb= taob.*(r-1).^2./(1-r.*taob.*p).*2 ;
Tdc= taoc.*(r-1).^2./(1-r.*taoc.*p).*2 ;
%ttt = (abs((1-r).^2.*taob.^0.5.*p./(1-r.*taob.*p))).^2;
plot(v,Tda/(1.096e-7),'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,Tdb/(1.096e-7),'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,Tdc/(1.096e-7),'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
legend('pump=xx','pump=x','pump=xxx')
%plot(v,TTT,'g','linewidth',2); xlabel('Frequency(Hz）');ylabel('Tran_Response');hold on;
%%理想透射谱
% n1 = 0.3;
% n2 = 0.6;
% n3 = 1;
% %H1_iii = abs(1./(1i*2*pi*v+0.5)).^0.6;
% H1_ideal=1./(abs(1i*2*pi*v*1e-9)).^(2*n1);
% H2_ideal=(1./(abs(1i*2*pi*v*1e-9)).^(2*n2));
% H3_ideal=(1./(abs(1i*2*pi*v*1e-9)).^(2*n3));
% figure;
% subplot(1,3,1);
% plot(v,H1_ideal,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
% subplot(1,3,2);
% plot(v,H2_ideal,'g','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H2');hold on;
% subplot(1,3,3);
% plot(v,H3_ideal,'b','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H3');hold on;
figure;
kk1 = 1./(abs(2*pi*(v+0.4e9)*1e-9)+0.003).^1;
kk2 = 1./(abs(2*pi*v*0.6e-9)+0.004).^1;
kk3 = 1./(abs(2*pi*v*1e-9)+0.006).^1;
plot(v,kk1/1.248,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
plot(v,kk2/1.248,'g','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
plot(v,kk3/1.248,'b','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission H1');hold on;
legend('pump=0.33','pump=.48','pump=.57')
%可以通过改变v的取值范围达到同样的效果，以上

%%理想微分方程结果




%% 输出结果
figure;
H_res1 = gauss_spec .* kk1;

idea_gauss_diff_power_spec=(abs(H_res1)).^2;
subplot(1,2,1)
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');title('理想微分功率谱**VI**');hold on;
subplot(1,2,2)
%%%微分输出结果
plot(double_f*1e-9,abs(H_res1)/0.02,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('idea_gauss_diff(a.u.)');hold on;

idea_t_c = fftshift(ifftshift(H_res1));
figure;
plot(t*1e+9,idea_t_c,'linewidth',2.5);



%% 微环输出
figure;
Tdaa= taoa.*(r-1)./(1-r.*taoa.*p) ;
H_ring_res1 = gauss_spec .* Tdaa;
ring_gauss_diff_power_spec=(abs(H_ring_res1)).^2;
subplot(1,2,1)
plot(double_f*1e-9,ring_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');title('理想微分功率谱**VI**');hold on;
subplot(1,2,2)
%%%微分输出结果
plot(double_f*1e-9,abs(H_ring_res1),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('idea_gauss_diff(a.u.)');hold on;

idea_t_cc = fftshift(ifftshift(H_ring_res1));
figure;
plot(t*1e+9,abs(idea_t_cc),'linewidth',2.5);

