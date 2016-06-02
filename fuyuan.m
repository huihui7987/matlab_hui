clc;
clear;
FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=601;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

n1=1.5;
n2=0.66;
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
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');title('����΢�ֹ�����**VI**');hold on;
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
plot(double_f*1e-9,Hr1,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('ģ��΢�����亯��');hold on;
subplot(1,3,2)
plot(double_f*1e-9,Hr1_Pi,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('ģ��΢����λ��');hold on;
subplot(1,3,3)
plot(double_f*1e-9,Hr1_tran,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('ģ��΢��͸����');hold on;


%%%%%%%��˹��������ģ��΢��΢����֮��%%%%%%%%%%
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
plot(double_f*1e-9,MRR_gauss_diff_power_spec,'g','linewidth',2.5);title('ģ��΢�����书����');hold on;
subplot(1,2,2)
plot(double_f*1e-9,MRR_gauss_diff,'g','linewidth',2.5);title('ģ��΢�����Ƶ��');hold on;


%%%%%%%%%%��˹��������ʵ��΢��֮��%%%%%%%%�˴����������⣬��ʱ��û�б�¶20160531_23:10

R=50e-6;
lamda=1439.2e-9:1e-12:1439.8e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
neff=3.17995709;
%neff=3.17996;
r=0.83;

L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
taoa=0.83;
Ta1= exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));
%Ta1= (taoa-r)./(1-r.*taoa);
T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r<=taoa
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end
%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))
taob=0.801;
Ta2= exp(1i*(pi+phi)).*(taob-r.*exp(-1i*phi))./(1-r.*taob.*exp(1i*phi));
T2= (abs(Ta2)).^2;%~~~~-
PHI2 = angle(Ta2);%~~

taoc=0.846;%1.5�ף���Ӧn1=1.5
Ta3= exp(1i*(pi+phi)).*(taoc-r.*exp(-1i*phi))./(1-r.*taoc.*exp(1i*phi));
T3= (abs(Ta3)).^2;%~~~~-
PHI3 = angle(Ta3);%~~
if r<=taoa
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end

ffff =Ta3 .* gauss_spec;

figure;
subplot(1,2,1);
plot(double_f,T1,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
legend('Critical-coupled','Under-coupled','Over-coupled')
subplot(1,2,2);
plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
legend('Critical-coupled','Under-coupled','Over-coupled')
%subplot(1,3,3);plot(v,y3,'linewidth',2);title('n=1.5');
figure;
%subplot(1,3,3);
plot(double_f,abs(ffff).^2/4/1e-6,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('res');hold on;

figure;
subplot(1,3,1);
plot(v,Ta3,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Ta3');hold on;
subplot(1,3,2);
plot(v,gauss_spec,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('gauss_spec');hold on;

figure;
H2=(1i*2*pi*double_f*1e-9).^(n2);
%����΢�ֽ��
idea_gauss_diff_v=gauss_spec.*H2;
%����΢�ֹ�����
idea_gauss_diff_power_spec_v=(abs(idea_gauss_diff_v)).^2;

plot(double_f*1e-9,idea_gauss_diff_power_spec/4.84,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('����΢��1.5�׹�����');
plot(double_f*1e-9,MRR_gauss_diff_power_spec/4.86,'g','linewidth',2.5);xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('ģ��΢�����书����');
plot(double_f*1e-9,abs(ffff).^2/4/1e-6/0.996,'b','linewidth',2.5);xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('MMR_based΢�����������');
plot(double_f*1e-9,idea_gauss_diff_power_spec_v/7.132/1e-3,'y','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('����΢��0.66�׹�����');
plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');%title('�����˹���幦����');
legend('ideal n=1.5','mod n=1.5','Ring based n=1.5','Ring based n=0.66','Input')












