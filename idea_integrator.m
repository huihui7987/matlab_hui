clear all
% lamda=1562.7e-9:1e-12:1563.1e-9;
% v=(3e8./lamda)-(3e8./1562.9e-9);

lamda=1538.7e-9:1e-12:1540.3e-9;
v=(3e8./lamda)-(3e8./1539.5e-9);

%%%%%%%%%
%����
%%%%%%%%%
n1=0.3;
n2=0.6;
n3=0.9;
%p1=exp(-1i*pi/2);
%p2=exp(1i*pi/2);
%t=2*pi*v;

%%%%%%%%%%����͸����%%%%%%%%%%%%
H1=(1./(abs(1i*2*pi*v)).^(2*n1));
H2=(1./(abs(1i*2*pi*v)).^(2*n2));
H3=(1./(abs(1i*2*pi*v)).^(2*n3));

h1=mapminmax(abs(H1),0,1);
h2=mapminmax(abs(H2),0,1);
h3=mapminmax(abs(H3),0,1);

subplot(1,3,1);
plot(v,H1,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
subplot(1,3,2);
plot(v,H2,'g','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H2');hold on;
subplot(1,3,3);
plot(v,H3,'b','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H3');hold on;
%legend('0.5','1','1.5')
% figure;
% %%%%%%%%������λ��Ӧ%%%%%%%%%%%%%%
% p1=(n1*pi/2).*(v>0)+(-n1*pi/2).*(v<0);
% p2=(n2*pi/2).*(v>0)+(-n2*pi/2).*(v<0);
% p3=(n3*pi/2).*(v>0)+(-n3*pi/2).*(v<0);
% 
% plot(v,p1,'r','linewidth',2.5);hold on;
% plot(v,p2,'g','linewidth',2.5); hold on;
% plot(v,p3,'b','linewidth',2.5);hold on;
% % plot(v,p1,'r');hold on;
% % plot(v,p2,'g'); hold on;
% % plot(v,p3,'b');hold on;

figure;
%%%%%%%%%%%%%%%%���봫�亯��%%%%%%%%%%%%%%%%%

Hv1=(1./(1i*2*pi*v)).^(n1);
Hv2=(1./(1i*2*pi*v)).^(n2);
Hv3=(1./(1i*2*pi*v)).^(n3);
subplot(1,3,1)
plot(v,Hv1,'r'); xlabel('Frequency(Hz��');ylabel('yyyyy');
subplot(1,3,2)
plot(v,Hv2,'r'); xlabel('Frequency(Hz��');ylabel('yyyyy');
subplot(1,3,3)
plot(v,Hv3,'r'); xlabel('Frequency(Hz��');ylabel('yyyyy');
figure;
pp1=angle(Hv1);
pp2=angle(Hv2);
pp3=angle(Hv3);
subplot(1,3,1)
plot(v,pp1,'r','linewidth',2.5); xlabel('xxxxx');ylabel('xxxxxxxxxxxxxxxx');hold on;
subplot(1,3,2)
plot(v,pp2,'r','linewidth',2.5); xlabel('xxxxx');ylabel('xxxxxxxxxxxxxxxx');hold on;
subplot(1,3,3)
plot(v,pp3,'r','linewidth',2.5); xlabel('xxxxx');ylabel('xxxxxxxxxxxxxxxx');hold on;

% figure;
% 
% H1=(abs(2*pi*v)).*(n1)*exp(1i*n1*pi/2).*(v>0)+(abs(2*pi*v)).*(n1)*exp(-1i*n1*pi/2).*(v<0);
% H2=(abs(2*pi*v)).*(n2)*exp(1i*n2*pi/2).*(v>0)+(abs(2*pi*v)).*(n2)*exp(-1i*n2*pi/2).*(v<0);
% H3=(abs(2*pi*v)).*(n3)*exp(1i*n3*pi/2).*(v>0)+(abs(2*pi*v)).*(n3)*exp(-1i*n3*pi/2).*(v<0);
% subplot(1,3,1)
% plot(v,H1,'r'); xlabel('Frequency(Hz��');ylabel('xxx');hold on;
% subplot(1,3,2)
% plot(v,H2,'g'); xlabel('Frequency(Hz��');ylabel('xxx');hold on;
% subplot(1,3,3)
% plot(v,H3,'b'); xlabel('Frequency(Hz��');ylabel('xxx');hold on;
% legend('0.5','1','1.5')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ģ��Ӧ�����ⲻ��
figure;
t0=5e-13;
Hw = exp(1i*v*t0).*Hv1;
HHH = angle(Hw);
subplot(1,2,1);
plot(v,abs(Hw),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
subplot(1,2,2);
plot(v,HHH,'g','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%5%%%%%
FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=1601;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

% n1=1.5;
% n2=0.66;
% n3=1.5;

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
%����΢�ֽ��
idea_gauss_diff=gauss_spec.*Hv1;
%����΢�ֹ�����
idea_gauss_diff_power_spec=(abs(idea_gauss_diff)).^2;
subplot(1,2,1)
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');title('����΢�ֹ�����**VI**');hold on;
subplot(1,2,2)
%%%΢��������
plot(double_f*1e-9,abs(idea_gauss_diff),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('idea_gauss_diff(a.u.)');hold on;









