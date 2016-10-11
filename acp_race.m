clc;
clear;
FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=1601;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

R=50e-6;
%lamda=1545.45e-9:1e-12:1546.15e-9;
lamda=1545e-9:1e-12:1546.6e-9;
v=(3e8./lamda)-(3e8./1545.8e-9);
%neff=3.1799883;
neff = 3.1799882;%小，右边高
r=0.88;
%yt=0.999;
Lc = R;
L = 2*pi*R+2*Lc;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi);
taoa=0.88;
yt=0.999;

% n1=1.5;
% n2=0.66;
% n3=1.5;

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');

%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
%gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行移位操作。
gauss_spec=fftshift(fft((gauss_time)));
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
%传输函数,实现n1 = 1.5阶
no=1.5;%此处可以修改n1的值，得到不同阶数的理想微分输出结果
nu=0.4;
nc=1;
H_idea=(1i*2*pi*double_f*1e-9).^(no);

%%%分段，暂时不知道错误在哪里，目前确定不对；0512
%H1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
idea_Intensity_Transmission = mapminmax((abs(1i*2*pi*double_f*1e-9)).^(2*no),0,1);
% idea_power_spec = (abs(H1)).^2;
 idea_phase_spec = angle(H_idea);
 
%%%%%作图%%%%%%
%%%%%%%%%%%%%%%
 %理想传输函数
 subplot(1,3,1)
 plot(double_f*1e-9,abs(H_idea),'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('idea_tran_spec');title('n1阶理想传输函数');hold on;
 %理想传输函数相位
 subplot(1,3,2)
 plot(double_f*1e-9,idea_phase_spec,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('idea_phase_spec');title('理想相位变化');hold on;
 %理想传输函数透射谱
subplot(1,3,3)
plot(double_f*1e-9,idea_Intensity_Transmission,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');title('理想透射谱');hold on;


%%%%%高斯函数理想微分之后的功率谱%%%%%
figure;
%理想微分结果
idea_gauss_diff=gauss_spec.*H_idea;
%理想微分功率谱
idea_gauss_diff_power_spec=(abs(idea_gauss_diff)).^2;
subplot(1,2,1)
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');title('理想微分功率谱**VI**');hold on;
subplot(1,2,2)
%%%微分输出结果
plot(double_f*1e-9,abs(idea_gauss_diff),'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('idea_gauss_diff(a.u.)');hold on;


%%微环微分器模型,c,o,u,分别代表严格耦合，过耦合和欠耦合

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%以下实现1.5阶%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=-3*1e-4;
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
Ho=(-1i*2*pi*double_f*1e-9).^(no);
Hro= -exp(-1i*double_f*1e-9*t0).*Ho;
Hro_tran = abs(Hro).^2;
Hro_Pi = angle(Hro);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%以上实现1.5阶%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%以下实现0.4阶%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% t0=-3*1e-4;
% n9=0.4;
% Ho=(-1i*2*pi*double_f*1e-9).^(n9);
% Hro= exp(-1i*double_f*1e-9*t0).*Ho;
% Hro_tran = abs(Hro).^2;
% Hro_Pi = angle(Hro);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%以上实现0.4阶%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%微环模型传输函数
figure;
subplot(1,3,1)
plot(double_f*1e-9,Hro,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('MRR_tran_spec');title('模型微环传输函数');hold on;
subplot(1,3,2)
plot(double_f*1e-9,Hro_Pi,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('MRR_tran_spec');title('模型微环相位谱');hold on;
subplot(1,3,3)
plot(double_f*1e-9,Hro_tran,'r','linewidth',2.5); xlabel('Frequency(Hz）');ylabel('MRR_tran_spec');title('模型微环透射谱');hold on;


%%%%%%%高斯函数经过模型微环微分器之后%%%%%%%%%%
figure;
%t0=-3*1e-4;
%H1=(1i*2*pi*double_f*1e-9).^(n1);
%H3=(-1i*2*pi*double_f*1e-9).^(n1);
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
%传输函数
%Hr1= -exp(-1i*double_f*1e-9*t0).*H3;
%微分结果
MRR_gauss_diff = gauss_spec.*Hro;
%微环微分之后功率谱
MRR_gauss_diff_power_spec=abs(MRR_gauss_diff).^2;

subplot(1,2,1)
plot(double_f*1e-9,MRR_gauss_diff_power_spec,'g','linewidth',2.5);title('模型微环传输功率谱');hold on;
subplot(1,2,2)
plot(double_f*1e-9,abs(MRR_gauss_diff),'g','linewidth',2.5);title('模型微环输出频谱');hold on;


%%%%%%%%%%高斯函数经过实际微环之后%%%%%%%%此处或许还有问题，暂时还没有暴露20160531_23:10

% R=50e-6;
% lamda=1439.2e-9:1e-12:1439.8e-9;
% v=(3e8./lamda)-(3e8./1439.5e-9);
% %neff=3.179992;
% neff = 3.15;
% r=0.98;
% %yt=0.999;
% Lc = R;
% L = 2*pi*R+2*Lc;
% phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
% p=exp(1i*phi);
% taoa=0.98;
% yt=0.999;

% R=50e-6;
% %lamda=1545.45e-9:1e-12:1546.15e-9;
% % lamda=1545.2e-9:1e-12:1546.4e-9;
% v=(3e8./lamda)-(3e8./1545.8e-9);
% neff=3.1799883;
% %neff = 3.17999;
% r=0.98;
% %yt=0.999;
% Lc = R;
% L = 2*pi*R+2*Lc;
% phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
% p=exp(1i*phi);
% taoa=0.98;
% yt=0.999;

%Ta1= exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));
Ta1= (r-taoa.*p)./(1-r.*taoa.*p);%%%%%%%%
%Ta1 = r*((1-taoa*yt^4.*p)./(1-r^2*taoa*yt^4.*p));
%Ta1= exp(1i*(pi+phi)).*(taoa*(yt)-r.*exp(-1i*phi))./(1-r.*taoa*(yt).*exp(1i*phi));

%Ta1=r*((1-yt*p)/(1-r^2*yt*p))


T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r<=taoa
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end
%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))

taob=0.875;
Ta2= (r-taob.*p)./(1-r.*taob.*p);
%Ta2 = r*((1-taob*yt^4.*p)./(1-r^2*taob*yt^4.*p));
T2= (abs(Ta2)).^2;
PHI2 = angle(Ta2);
if r<=taob
    PHI2 = PHI2+(PHI2<0)*2*pi ;
end

taoc=0.89;
Ta3= (r-taoc.*p)./(1-r.*taoc.*p);
%Ta3 = r*((1-taoc*yt^4.*p)./(1-r^2*taoc*yt^4.*p));
T3= (abs(Ta3)).^2;
PHI3 = angle(Ta3);
if r<=taoc
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end

ff_o =Ta3 .* gauss_spec;
ff_u = Ta2 .* gauss_spec;

figure;
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
plot(double_f,abs(ff_o).^2/4/1e-6,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('res');hold on;

figure;
subplot(1,3,1);
plot(v,Ta3,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Ta3');hold on;
subplot(1,3,2);
plot(v,gauss_spec,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('gauss_spec');hold on;

figure;
n2 = 0.66;
H2=(1i*2*pi*double_f*1e-9).^(n2);
%理想微分结果
idea_gauss_diff_v=gauss_spec.*H2;
%理想微分功率谱
idea_gauss_diff_power_spec_v=(abs(idea_gauss_diff_v)).^2;

plot(double_f*1e-9,idea_gauss_diff_power_spec/4.84,'r','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('理想微分1.5阶功率谱');
%plot(double_f*1e-9,MRR_gauss_diff_power_spec/4.86,'g','linewidth',2.5);xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('模型微环传输功率谱');
plot(double_f*1e-9,abs(ff_o).^2/4/1e-6/0.996/18.84,'b','linewidth',2.5);xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('MMR_based微分输出功率谱');
%plot(double_f*1e-9,idea_gauss_diff_power_spec_v/7.132/1e-3,'y','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('理想微分0.66阶功率谱');
plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');%title('输入高斯脉冲功率谱');
legend('ideal n=1.5','mod n=1.5','Ring based n=1.5','Ring based n=0.66','Input')

figure;
%%%%微环一阶微分与输入脉冲同图
H_ide_1st = 1i*2*pi*double_f*1e-9;
H_1st = -exp(-1i*v*t0).*H_ide_1st;%建模型
H_1st_idel = gauss_spec.*H_ide_1st;%理想输出
H_1st_Mrr = gauss_spec.*Ta1;%Mrr输出
H_1st_res = gauss_spec.*H_1st;%模型输出
plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('输入高斯脉冲功率谱');
plot(double_f*1e-9,abs(H_1st_idel).^2/0.09218,'r','linewidth',2.5);xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('理想一阶微分输出功率谱');
plot(double_f*1e-9,abs(H_1st_Mrr).^2/(6.325e-5),'g','linewidth',2.5);xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('MMR一阶微分输出功率谱');
%plot(double_f*1e-9,abs(H_1st_res).^2/0.09243,'b','linewidth',2.5);xlabel('Frequency(GHz）');ylabel('power(a.u.)');hold on;%title('模型一阶微分输出功率谱');
legend('Input','ideal n=1','Ring based n=1');

%%%%%%%%%%%%%%%%%%%%时域验证
figure;
Mrr_t_out =ifft(H_1st_Mrr);

ide_t_out = ifft(H_1st_idel);%1阶理想
plot(t*1e+9,abs(Mrr_t_out),'r','linewidth',2); xlabel('Time/ns');ylabel('Amplitude/V');title('微环微分时域');hold on;
figure;
plot(t*1e+9,abs(ide_t_out),'g','linewidth',2); xlabel('Time/ns');ylabel('Amplitude/V');title('理想微分时域');hold on;

%%%%%%%%%%%%%%%%%相位验证
figure;

plot(v,PHI2,'r','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;
phi_under_out = angle((-1i*2*pi*double_f*1e-9).^(0.7));
plot(v,phi_under_out,'b','linewidth',2); xlabel('Frequency(Hz）');ylabel('Phase Response');hold on;















