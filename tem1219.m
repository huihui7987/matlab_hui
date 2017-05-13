clc;
clear;
FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=601;                %������
dt=time_window/(Ns-1);  %����ʱ����
t=0:dt:time_window;     %����ʱ��

kk = fftshift((2*pi/time_window).*[(-Ns/2):(Ns/2-1)]);
lambda0 = 1450e-9;
f0=3e8./lambda0;
omega0 = 2*pi*f0;
lamd = 2*pi*3e8./(kk+omega0);


% n1=1.5;
% n2=0.66;
% n3=1.5;

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');

%===========���¼���˫���ס�˫�߹����ס�˫�߹������ܶ�=================
%gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %����Ҷ�任�����ҽ�����λ������
gauss_spec=fftshift(ifft(gauss_time)); 
gauss_spec=gauss_spec/Ns;   %��ʵ�ʵķ���ֵ����һ����
df=1/time_window;               %Ƶ�ʷֱ���
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %˫��Ƶ�׶�Ӧ��Ƶ��

%figure;

%plot(double_f*1e-9,gauss_spec,'linewidth',2.5);
figure; %������
plot(fftshift(lamd),abs(real(gauss_spec)),'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Amplitude/V');
title('double Amplitude spectrum');
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
plot(fftshift(lamd),double_power_spec_mW/0.1132,'linewidth',2.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');
figure;
%plot(fftshift(lamd),Aw_dB,'b--','linewidth',2);hold on;

%%%%%%%%%%%%%%%%%
%%����΢�������亯��
%%%ʵ��֤���ú������ã�Ϊ������ʽ��matlab��ͼֻ�ܻ�����ʵ������������λ�ɷֱ�����abs()��angle()�鿴��0512
%%%%%%%%%%%%%%%%%
%���亯��,ʵ��n1 = 1.5��
no=1.5;%�˴������޸�n1��ֵ���õ���ͬ����������΢��������
nu=0.4;
nc=1;
H_idea=(1i*2*pi*double_f*1e-9).^(no);

%%%�ֶΣ���ʱ��֪�����������Ŀǰȷ�����ԣ�0512
%H1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
idea_Intensity_Transmission = mapminmax((abs(1i*2*pi*double_f*1e-9)).^(2*no),0,1);
% idea_power_spec = (abs(H1)).^2;
 idea_phase_spec = angle(H_idea);
 
%%%%%��ͼ%%%%%%
%%%%%%%%%%%%%%%
 %���봫�亯��
 subplot(1,3,1)
 plot(double_f*1e-9,abs(H_idea),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('idea_tran_spec');title('n1�����봫�亯��');hold on;
 %���봫�亯����λ
 subplot(1,3,2)
 plot(double_f*1e-9,idea_phase_spec,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('idea_phase_spec');title('������λ�仯');hold on;
 %���봫�亯��͸����
subplot(1,3,3)
plot(double_f*1e-9,idea_Intensity_Transmission,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');title('����͸����');hold on;


%%%%%��˹��������΢��֮��Ĺ�����%%%%%
figure;
%����΢�ֽ��
idea_gauss_diff=gauss_spec.*H_idea;
%����΢�ֹ�����
idea_gauss_diff_power_spec=(abs(idea_gauss_diff)).^2;
subplot(1,2,1)
plot(double_f*1e-9,idea_gauss_diff_power_spec,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');title('����΢�ֹ�����**VI**');hold on;
subplot(1,2,2)
%%%΢��������
plot(double_f*1e-9,abs(idea_gauss_diff),'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('idea_gauss_diff(a.u.)');hold on;


%%΢��΢����ģ��,c,o,u,�ֱ�����ϸ���ϣ�����Ϻ�Ƿ���

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%����ʵ��1.5��%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=-3*1e-4;
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
Ho=(-1i*2*pi*double_f*1e-9).^(no);
Hro= -exp(-1i*double_f*1e-9*t0).*Ho;
Hro_tran = abs(Hro).^2;
Hro_Pi = angle(Hro);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%����ʵ��1.5��%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%����ʵ��0.4��%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% t0=-3*1e-4;
% n9=0.4;
% Ho=(-1i*2*pi*double_f*1e-9).^(n9);
% Hro= exp(-1i*double_f*1e-9*t0).*Ho;
% Hro_tran = abs(Hro).^2;
% Hro_Pi = angle(Hro);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%����ʵ��0.4��%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%΢��ģ�ʹ��亯��
figure;
subplot(1,3,1)
plot(double_f*1e-9,Hro,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('ģ��΢�����亯��');hold on;
subplot(1,3,2)
plot(double_f*1e-9,Hro_Pi,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('ģ��΢����λ��');hold on;
subplot(1,3,3)
plot(double_f*1e-9,Hro_tran,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('MRR_tran_spec');title('ģ��΢��͸����');hold on;


%%%%%%%��˹��������ģ��΢��΢����֮��%%%%%%%%%%
figure;
%t0=-3*1e-4;
%H1=(1i*2*pi*double_f*1e-9).^(n1);
%H3=(-1i*2*pi*double_f*1e-9).^(n1);
%H3=(abs(-2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(-2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
%���亯��
%Hr1= -exp(-1i*double_f*1e-9*t0).*H3;
%΢�ֽ��
MRR_gauss_diff = gauss_spec.*Hro;
%΢��΢��֮������
MRR_gauss_diff_power_spec=abs(MRR_gauss_diff).^2;

subplot(1,2,1)
plot(double_f*1e-9,MRR_gauss_diff_power_spec,'g','linewidth',2.5);title('ģ��΢�����书����');hold on;
subplot(1,2,2)
plot(double_f*1e-9,MRR_gauss_diff,'g','linewidth',2.5);title('ģ��΢�����Ƶ��');hold on;


%%%%%%%%%%��˹��������ʵ��΢��֮��%%%%%%%%�˴����������⣬��ʱ��û�б�¶20160531_23:10

R=50e-6;
lamda=1439.2e-9:1e-12:1439.8e-9;
v=(3e8./lamda)%-(3e8./1439.5e-9);
neff=3.17995709;
%neff=3.17996;
r=0.83;

L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
p=exp(1i*phi/2);

taoa=0.83;
%Ta1=
%exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));%����Դ��ʽ������
%Ta1= (taoa-r)./(1-r.*taoa);
Ta1 = (r-taoa*p.^2)./(1-taoa*r*p.^2);%���ƹ�ʽ������
T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r<=taoa
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end
%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))
taob=0.801;%��Ӧ0.4��
%Ta2= exp(1i*(pi+phi)).*(taob-r.*exp(-1i*phi))./(1-r.*taob.*exp(1i*phi));
Ta2 = (r-taob*p.^2)./(1-taob*r*p.^2);
T2= (abs(Ta2)).^2;%~~~~-
PHI2 = angle(Ta2);%~~

taoc=0.846;%1.5�ף���Ӧn1=1.5
%Ta3= exp(1i*(pi+phi)).*(taoc-r.*exp(-1i*phi))./(1-r.*taoc.*exp(1i*phi));
Ta3 = (r-taoc*p.^2)./(1-taoc*r*p.^2);
T3= (abs(Ta3)).^2;%~~~~-
PHI3 = angle(Ta3);%~~
if r<=taoa
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end

ff_o =Ta3 .* gauss_spec;
ff_u = Ta2 .* gauss_spec;
ff_c = Ta1 .*gauss_spec;
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
plot(double_f,abs(ff_o).^2/4/1e-6,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('res');hold on;

figure;
subplot(1,3,1);
plot(v,Ta3,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Ta3');hold on;
subplot(1,3,2);
plot(v,gauss_spec,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('gauss_spec');hold on;

figure;
n2 = 0.66;
H2=(1i*2*pi*double_f*1e-9).^(n2);
%����΢�ֽ��
idea_gauss_diff_v=gauss_spec.*H2;
%����΢�ֹ�����
idea_gauss_diff_power_spec_v=(abs(idea_gauss_diff_v)).^2;
plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('�����˹���幦����');
plot(double_f*1e-9,idea_gauss_diff_power_spec/4.84,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('����΢��1.5�׹�����');
%plot(double_f*1e-9,MRR_gauss_diff_power_spec/4.86,'g','linewidth',2.5);xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('ģ��΢�����书����');
plot(double_f*1e-9,abs(ff_o).^2/4/1e-6/0.996,'b','linewidth',2.5);xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('MMR_based΢�����������');
%plot(double_f*1e-9,idea_gauss_diff_power_spec_v/7.132/1e-3,'y','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('����΢��0.66�׹�����');
%plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');%title('�����˹���幦����');
legend('Input','ideal n=1.5','Ring based n=1.5')
%%%%%%%%%%%%%%%%%%%%%%1��
figure;
nn=1;
H_idea_c=(1i*2*pi*double_f*1e-9).^(nn);%1������
idea_gauss_diff_c=gauss_spec.*H_idea_c;
%i������΢�ֹ�����
idea_gauss_diff_power_spec_c=(abs(idea_gauss_diff_c)).^2;
%1��΢��

plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('�����˹���幦����');
plot(double_f*1e-9,idea_gauss_diff_power_spec_c/0.09242,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('����΢��1�׹�����');
plot(double_f*1e-9,abs(ff_c).^2/(3.506e-6),'b','linewidth',2.5);xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('MMR_based΢�����������');
legend('Input','ideal n=1','Ring based n=1')

%%%%%%%%%%%%%0.4��
figure;
nm=0.4;
H_idea_u=(1i*2*pi*double_f*1e-9).^(nm);%0.4������
idea_gauss_diff_u=gauss_spec.*H_idea_u;
%i������΢�ֹ�����
idea_gauss_diff_power_spec_u=(abs(idea_gauss_diff_u)).^2;
%1��΢��

plot(double_f*1e-9,double_power_spec_mW/0.1132,'k','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('�����˹���幦����');
plot(double_f*1e-9,idea_gauss_diff_power_spec_u/0.001145,'r','linewidth',2.5); xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('����΢��1�׹�����');
plot(double_f*1e-9,abs(ff_u).^2/(3.281e-6),'b','linewidth',2.5);xlabel('Frequency(GHz��');ylabel('power(a.u.)');hold on;%title('MMR_based΢�����������');
legend('Input','ideal n=0.4','Ring based n=0.4')

%%%%%%%%%%%%%ʱ����֤%%%%%%%%
%%%1��%%%%%%%%
%%1������
figure;
idea_t_c = ifftshift(idea_gauss_diff_c);
plot(t*1e+9,idea_t_c,'linewidth',2.5);





