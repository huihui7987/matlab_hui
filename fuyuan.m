clc;
clear;
FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=2048;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

n1=1.5;
n2=1;
n3=1.5;

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');

%===========���¼���˫���ס�˫�߹����ס�˫�߹������ܶ�=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %����Ҷ�任�����ҽ�����λ������
gauss_spec=gauss_spec/Ns;   %��ʵ�ʵķ���ֵ����һ����
df=1/time_window;               %Ƶ�ʷֱ���
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %˫��Ƶ�׶�Ӧ��Ƶ��

% figure;
% 
% plot(double_f*1e-9,gauss_spec,'linewidth',2.5);
% figure; %������
% plot(double_f*1e-9,abs(gauss_spec),'linewidth',2.5);
% xlabel('Frequency/GHz');
% ylabel('Amplitude/V');
% title('double Amplitude spectrum');
% 
% 
% figure; %��λ��
% plot(double_f*1e-9,angle(gauss_spec),'linewidth',2.5);
% xlabel('Frequency/GHz');
% ylabel('Phase/rad');
% title('double Phase spectrum');

figure; %��˹���幦����
double_power_spec_W=abs(gauss_spec).^2;                 %˫�߹����ף���λW��
double_power_spec_mW=double_power_spec_W*1e+3;          %˫�߹����ף���λmW��
double_power_spec_dBm=10*log10(double_power_spec_mW);   %˫�߹����ף���λdBm��
plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');
figure;


%%%%%%%%%%%%%%%%%
%%����΢�������亯��
%%%ʵ��֤���ú������ã�Ϊ������ʽ��matlab��ͼֻ�ܻ�����ʵ������������λ�ɷֱ�����abs()��angle()�鿴��0512
%%%%%%%%%%%%%%%%%
%���亯��
H1=(1i*2*pi*double_f*1e-9).^(n1);
%%%�ֶΣ���ʱ��֪�����������Ŀǰȷ�����ԣ�0512
%H1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
idea_Intensity_Transmission = mapminmax((abs(1i*2*pi*double_f*1e-9)).^(2*n1),0,1);
% idea_power_spec = (abs(H1)).^2;
 idea_phase_spec = angle(H1);
 
%%%%%��ͼ%%%%%%
%%%%%%%%%%%%%%%
 %���봫�亯��
 subplot(1,3,1)
 plot(double_f*1e-9,abs(H1),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('idea_tran_spec');title('���봫�亯��');hold on;
 %���봫�亯����λ
 subplot(1,3,2)
 plot(double_f*1e-9,idea_phase_spec,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('idea_phase_spec');title('������λ�仯');hold on;
 %���봫�亯��͸����
subplot(1,3,3)
plot(double_f*1e-9,idea_Intensity_Transmission,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');title('����͸����');hold on;


%%%%%��˹��������΢��֮��Ĺ�����%%%%%
figure;
%����΢�ֽ��
idea_gauss_diff=gauss_spec.*H1;
%����΢�ֹ�����
idea_gauss_diff_power_spec=(abs(idea_gauss_diff)).^2;
subplot(1,2,1)
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');title('����΢�ֹ�����');hold on;
subplot(1,2,2)
%%%΢��������
plot(double_f*1e-9,abs(idea_gauss_diff),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('idea_gauss_diff(a.u.)');hold on;


%%΢��΢�������亯��
%n=1.5
figure;
t0=-3*1e-4;
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
H3=(-1i*2*pi*double_f*1e-9).^(n1);
Hr1= -exp(-1i*double_f*1e-9*t0).*H3;
Hr1_tran = abs(Hr1).^2;
Hr1_Pi = angle(Hr1);

%΢�����亯��
subplot(1,3,1)
plot(double_f*1e-9,Hr1,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('΢�����亯��');hold on;
subplot(1,3,2)
plot(double_f*1e-9,Hr1_Pi,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');hold on;
subplot(1,3,3)
plot(double_f*1e-9,Hr1_tran,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');hold on;


%%%%%%%��˹��������΢��΢����֮��%%%%%%%%%%
figure;
%t0=-3*1e-4;
%H1=(1i*2*pi*double_f*1e-9).^(n1);
%H3=(-1i*2*pi*double_f*1e-9).^(n1);
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
%���亯��
%Hr1= -exp(-1i*double_f*1e-9*t0).*H3;
%΢�ֽ��
MRR_gauss_diff = gauss_spec.*Hr1;
%΢��΢��֮������
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
% plot(v,TT1,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Tran_Response');hold on;






