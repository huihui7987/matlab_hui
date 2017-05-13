clc;
clear;


FWHM=62e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=5000*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=201;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

% n1=1.5;
% n2=0.66;
% n3=1.5;
figure;
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

plot(double_f*1e-9,gauss_spec,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('Tdaa');title('gau');hold on;
%%

R=50e-6;
lamda=1439.4e-9:1e-12:1439.6e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
%neff=3.179962;
neff=3.40002;
r=0.9;
%yt=0.999;
Lc = R;
L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi/2);
taoa=0.99999;

figure;
%Tdaa= taoa.*(r-1)./(1-r.*taoa.*p) ;
r1 = 0.999;
r2 = 0.999;
k1 = sqrt(1-r1^2);
k2 = sqrt(1-r2^2);

Tdaa = (-taoa*k1*k2*p)./(1-taoa*r1*r2*p.^2);
H_ring_res1 = gauss_spec .* Tdaa;

plot(double_f*1e-9,abs(Tdaa),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('Tdaa');title('���亯��');hold on;

ring_gauss_diff_power_spec=(abs(H_ring_res1)).^2;
figure;
plot(double_f*1e-9,abs(H_ring_res1),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('Tdaa');title('Ƶ����');hold on;
figure;
hht = abs(ifft((H_ring_res1)))/1.886e-117/0.7889/1.733;
plot(t*1e+9,hht,'linewidth',2.5);


