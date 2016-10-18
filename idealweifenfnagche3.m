clc;
clear;


FWHM=800e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=401;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

% n1=1.5;
% n2=0.66;
% n3=1.5;
figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
plot(t*1e+9,gauss_time,'linewidth',2.5);
% %���ϣ���ͼ�����ۣ�ԭ��ȡ����̫��
% xlabel('Time/ns');
% ylabel('Amplitude/V');
% title('Gauss pulse');

%===========���¼���˫���ס�˫�߹����ס�˫�߹������ܶ�=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %����Ҷ�任�����ҽ�����λ������
%gauss_spec=fftshift(gauss_time); 
gauss_spec=gauss_spec/Ns;   %��ʵ�ʵķ���ֵ����һ����
df=1/time_window;               %Ƶ�ʷֱ���
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %˫��Ƶ�׶�Ӧ��Ƶ��
plot(double_f*1e-9,(gauss_spec),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('Tdaa');title('gauss_spec');hold

kk1 = 1./((1i*2*pi*double_f*1e-9)+0.048);
kk2 = 1./((1i*2*pi*double_f*1e-9)+0.028);
figure;
plot(double_f*1e-9,kk1,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('kk1-Intensity Transmission H1');hold on;
figure;
hh = gauss_spec .* kk1;
plot(double_f*1e-9,abs(hh),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
figure;
hht = ifft(hh);
plot(t*1e+9,abs(hht)/0.001953,'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');




