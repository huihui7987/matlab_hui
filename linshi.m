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
plot(v,PI,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
subplot(1,2,2)
plot(v,H1,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% plot(v,H2,'g'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% plot(v,H3,'b'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% legend('0.5','1','1.5')

figure;
H11=((abs(2*pi*v)).*(n1)).*exp(1i*n1*pi/2).*(v>0)+((abs(2*pi*v)).*(n1)).*exp(-1i*n1*pi/2).*(v<0);

plot(v,H11,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;



% 
% clc;
% clear;
% FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
% time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
% Ns=2048;                %������
% dt=time_window/(Ns-1);  %����ʱ����
% t=0:dt:time_window;     %����ʱ��
% 
% n1=1.5;
% n2=1;
% n3=1.5;
% 
% gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
% %plot(t*1e+9,gauss_time,'linewidth',2.5);
% xlabel('Time/ns');
% ylabel('Amplitude/V');
% title('Gauss pulse');
% 
% %===========���¼���˫���ס�˫�߹����ס�˫�߹������ܶ�=================
% gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %����Ҷ�任�����ҽ�����λ������
% gauss_spec=gauss_spec/Ns;   %��ʵ�ʵķ���ֵ��
% df=1/time_window;               %Ƶ�ʷֱ���
% k=floor(-(Ns-1)/2:(Ns-1)/2);    
% % k=0:Ns-1;
% double_f=k*df;   %˫��Ƶ�׶�Ӧ��Ƶ��
% 
% % figure;
% % 
% % plot(double_f*1e-9,gauss_spec,'linewidth',2.5);
% % figure; %������
% % plot(double_f*1e-9,abs(gauss_spec),'linewidth',2.5);
% % xlabel('Frequency/GHz');
% % ylabel('Amplitude/V');
% % title('double Amplitude spectrum');
% % 
% % 
% % figure; %��λ��
% % plot(double_f*1e-9,angle(gauss_spec),'linewidth',2.5);
% % xlabel('Frequency/GHz');
% % ylabel('Phase/rad');
% % title('double Phase spectrum');
% 
% %figure; %��˹���幦����
% double_power_spec_W=abs(gauss_spec).^2;                 %˫�߹����ף���λW��
% double_power_spec_mW=double_power_spec_W*1e+3;          %˫�߹����ף���λmW��
% double_power_spec_dBm=10*log10(double_power_spec_mW);   %˫�߹����ף���λdBm��
% %plot(double_f*1e-9,double_power_spec_mW/0.1132,'linewidth',2.5);
% xlabel('Frequency/GHz');
% ylabel('Power/dBm');
% title('double Power spectrum');
% 
% %%%%%%%%%%%%%%%%%
% %%����΢�������亯��
% %%%%%%%%%%%%%%%%%
% id_H1=(1i*2*pi*double_f*1e-9).^(n1);%%%��
% id_Hr1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
% subplot(1,2,1);
% plot(double_f*1e-9,id_H1,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% subplot(1,2,2)
% plot(double_f*1e-9,id_Hr1,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% 
% %%%΢����n=1.5
% 
% t0=-3*1e-4;
% Mr_H1=(exp(-1i*double_f*1e-9*t0)).*(-1i*2*pi*double_f*1e-9).^(n1);%%%��
% %%��
% H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
% Mr_Hr1= exp(-1i*double_f*1e-9*t0).*H3;
% 
% figure;
% subplot(1,2,1);
% plot(double_f*1e-9,Mr_H1,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% subplot(1,2,2);
% plot(double_f*1e-9,Mr_Hr1,'r'); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;






    
