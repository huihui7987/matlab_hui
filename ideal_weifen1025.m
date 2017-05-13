clc;
clear;


FWHM=62e-12;            %��˹�ź�FWHM��ȣ�Ϊ62ps
time_window=2000*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=1001;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

% n1=1.5;
% n2=0.66;
% n3=1.5;
% figure;
gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
% plot(t*1e+9,gauss_time,'linewidth',2.5);
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
%plot(double_f*1e-9,(gauss_spec),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('Tdaa');title('gauss_spec');hold

kk1 = 1./((1i*2*pi*double_f*1e-9)+0.038);
kk2 = 1./((1i*2*pi*double_f*1e-9)+0.054);
kk3 = 1./((1i*2*pi*double_f*1e-9)+0.070);
kk4 = 1./((1i*2*pi*double_f*1e-9)+0.084);
kk5 = 1./((1i*2*pi*double_f*1e-9)+0.11);
kk6 = 1./((1i*2*pi*double_f*1e-9)+0.13);
figure;
plot(double_f*1e-9,angle(kk1),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('kk1-Intensity Transmission H1');hold on;
figure;

res_rr = abs(kk6.^2);
plot(double_f*1e-9,res_rr/342.9,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('kk1-Intensity Transmission H1');hold on;

hh1 = gauss_spec .* kk1;
% plot(double_f*1e-9,abs(hh),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission H1');hold on;
figure;
hht1 = abs(ifft(hh1))/1.603e-6;
hh2 = gauss_spec .* kk2;
hht2 = abs(ifft(hh2))/1.603e-6;
hh3 = gauss_spec .* kk3;
hht3 = abs(ifft(hh3))/1.603e-6;
hh4 = gauss_spec .* kk4;
hht4 = abs(ifft(hh4))/1.603e-6;
hh5 = gauss_spec .* kk5;
hht5 = abs(ifft(hh5))/1.603e-6;
hh6 = gauss_spec .* kk6;
hht6 = abs(ifft(hh6))/1.603e-6;
plot(t*1e+9,abs(hht1),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht2),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht3),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht4),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht5),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
plot(t*1e+9,abs(hht6),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
legend('k=0.038','k=0.054','k=0.070','k=0.084','k=0.15','k=0.25')
% figure;
% hf1 = gauss_spec .* abs(kk1);
% hf2 = abs(ifft(hf1));
% hf3 = gauss_spec .* abs(kk5);
% hf4 = abs(ifft(hf3));
% plot(t*1e+9,abs(hf2),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
% plot(t*1e+9,abs(hf4),'linewidth',2.5);xlabel('Time(ps��');ylabel('Intensity(a.u.)');hold on;
